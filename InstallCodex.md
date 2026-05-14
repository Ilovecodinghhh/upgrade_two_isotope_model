# Installing Codex CLI on WSL with Mindflow API

> Codex CLI uses OpenAI's **Responses API** internally. Mindflow exposes a **Chat Completions API**. A lightweight Python proxy bridges the two.

---

## 1. Prerequisites

```bash
# Node.js 22+ (via nvm)
curl -o- https://raw.githubusercontent.com/nvm-sh/nvm/v0.40.3/install.sh | bash
source ~/.bashrc
nvm install 22
node -v  # should be v22.x+

# Python 3 + aiohttp (for the proxy)
sudo apt update && sudo apt install -y python3 python3-pip
pip3 install aiohttp
```

---

## 2. Install Codex CLI

```bash
npm install -g @openai/codex
codex --version  # should print 0.130.0+
```

---

## 3. Set Up the API Proxy

Codex speaks **Responses API** → proxy translates to **Chat Completions API** → Mindflow.

### 3a. Create the proxy script

```bash
mkdir -p ~/.codex
cat > ~/.codex/codex_proxy.py << 'PROXY_EOF'
#!/usr/bin/env python3
"""
codex_proxy.py — Translates OpenAI Responses API to Chat Completions API
for third-party providers (Mindflow, etc.)
"""
import asyncio, json, uuid, time, os, sys
from aiohttp import web, ClientSession, ClientTimeout

UPSTREAM_URL = os.environ.get("UPSTREAM_URL", "https://ai.mindflow.com.cn/v1/chat/completions")
LISTEN_HOST = "127.0.0.1"
LISTEN_PORT = 18900


def convert_input_to_messages(body):
    messages = []
    instructions = body.get("instructions", "")
    if instructions:
        messages.append({"role": "system", "content": instructions})

    inp = body.get("input", "")
    if isinstance(inp, str):
        messages.append({"role": "user", "content": inp})
    elif isinstance(inp, list):
        for item in inp:
            if isinstance(item, str):
                messages.append({"role": "user", "content": item})
            elif isinstance(item, dict):
                msg_type = item.get("type", "")
                if msg_type == "message":
                    role = item.get("role", "user")
                    content = item.get("content", "")
                    if isinstance(content, list):
                        parts = []
                        for part in content:
                            if isinstance(part, dict):
                                t = part.get("type", "")
                                if t in ("input_text", "output_text", "text"):
                                    parts.append(part.get("text", ""))
                            elif isinstance(part, str):
                                parts.append(part)
                        content = "\n".join(parts)
                    messages.append({"role": role, "content": content})
                elif msg_type == "function_call":
                    messages.append({
                        "role": "assistant",
                        "content": None,
                        "tool_calls": [{
                            "id": item.get("call_id", f"call_{uuid.uuid4().hex[:8]}"),
                            "type": "function",
                            "function": {
                                "name": item.get("name", ""),
                                "arguments": item.get("arguments", "{}"),
                            }
                        }]
                    })
                elif msg_type == "function_call_output":
                    messages.append({
                        "role": "tool",
                        "tool_call_id": item.get("call_id", ""),
                        "content": item.get("output", ""),
                    })
    return messages


def convert_tools(body):
    tools = []
    for tool in body.get("tools", []):
        if tool.get("type") == "function":
            tools.append({
                "type": "function",
                "function": {
                    "name": tool.get("name", ""),
                    "description": tool.get("description", ""),
                    "parameters": tool.get("parameters", {}),
                    "strict": tool.get("strict", False),
                }
            })
    return tools


async def call_upstream(api_key, body):
    model = body.get("model", "claude-opus-4-6")
    messages = convert_input_to_messages(body)
    tools = convert_tools(body)
    max_tokens = body.get("max_output_tokens", 16384)
    temperature = body.get("temperature", 1.0)

    chat_body = {
        "model": model,
        "messages": messages,
        "max_tokens": max_tokens,
        "temperature": temperature,
    }
    if tools:
        chat_body["tools"] = tools

    async with ClientSession() as session:
        async with session.post(
            UPSTREAM_URL,
            json=chat_body,
            headers={
                "Authorization": f"Bearer {api_key}",
                "Content-Type": "application/json",
            },
            timeout=ClientTimeout(total=300),
        ) as resp:
            return await resp.json()


def build_full_response(resp_id, chat_resp):
    choice = chat_resp.get("choices", [{}])[0]
    message = choice.get("message", {})
    usage = chat_resp.get("usage", {})
    output = []

    if message.get("tool_calls"):
        for tc in message["tool_calls"]:
            output.append({
                "type": "function_call",
                "id": f"fc_{uuid.uuid4().hex[:8]}",
                "call_id": tc.get("id", f"call_{uuid.uuid4().hex[:8]}"),
                "name": tc["function"]["name"],
                "arguments": tc["function"]["arguments"],
                "status": "completed",
            })
    else:
        output.append({
            "type": "message",
            "id": f"msg_{uuid.uuid4().hex[:8]}",
            "role": "assistant",
            "content": [{"type": "output_text", "text": message.get("content", "")}],
            "status": "completed",
        })

    return {
        "id": resp_id,
        "object": "response",
        "created_at": int(time.time()),
        "status": "completed",
        "model": chat_resp.get("model", "claude-opus-4-6"),
        "output": output,
        "usage": {
            "input_tokens": usage.get("prompt_tokens", 0),
            "output_tokens": usage.get("completion_tokens", 0),
            "total_tokens": usage.get("total_tokens", 0),
        },
    }


def build_stream_events(resp_id, chat_resp):
    events = []
    choice = chat_resp.get("choices", [{}])[0]
    message = choice.get("message", {})
    usage = chat_resp.get("usage", {})

    events.append({"type": "response.created", "response": {
        "id": resp_id, "object": "response", "status": "in_progress", "output": [],
    }})

    if message.get("tool_calls"):
        for i, tc in enumerate(message["tool_calls"]):
            item_id = f"fc_{uuid.uuid4().hex[:8]}"
            call_id = tc.get("id", f"call_{uuid.uuid4().hex[:8]}")
            fn = tc["function"]
            item = {
                "type": "function_call", "id": item_id, "call_id": call_id,
                "name": fn["name"], "arguments": "", "status": "in_progress",
            }
            events.append({"type": "response.output_item.added", "output_index": i, "item": item})
            events.append({"type": "response.function_call_arguments.delta",
                          "output_index": i, "item_id": item_id, "delta": fn["arguments"]})
            events.append({"type": "response.function_call_arguments.done",
                          "output_index": i, "item_id": item_id, "arguments": fn["arguments"]})
            item_done = dict(item)
            item_done["arguments"] = fn["arguments"]
            item_done["status"] = "completed"
            events.append({"type": "response.output_item.done", "output_index": i, "item": item_done})
    else:
        text = message.get("content", "")
        item_id = f"msg_{uuid.uuid4().hex[:8]}"
        events.append({"type": "response.output_item.added", "output_index": 0, "item": {
            "type": "message", "id": item_id, "role": "assistant", "content": [], "status": "in_progress",
        }})
        events.append({"type": "response.content_part.added", "output_index": 0, "content_index": 0,
                       "part": {"type": "output_text", "text": ""}})
        cs = 40
        for i in range(0, max(1, len(text)), cs):
            events.append({"type": "response.output_text.delta", "output_index": 0,
                          "content_index": 0, "delta": text[i:i+cs]})
        events.append({"type": "response.output_text.done", "output_index": 0,
                       "content_index": 0, "text": text})
        events.append({"type": "response.content_part.done", "output_index": 0, "content_index": 0,
                       "part": {"type": "output_text", "text": text}})
        events.append({"type": "response.output_item.done", "output_index": 0, "item": {
            "type": "message", "id": item_id, "role": "assistant",
            "content": [{"type": "output_text", "text": text}], "status": "completed",
        }})

    full = build_full_response(resp_id, chat_resp)
    events.append({"type": "response.completed", "response": full})
    return events


async def handle_responses(request):
    body = await request.json()
    api_key = request.headers.get("Authorization", "").replace("Bearer ", "")
    stream = body.get("stream", False)
    resp_id = f"resp_{uuid.uuid4().hex[:24]}"

    try:
        chat_resp = await call_upstream(api_key, body)
    except Exception as e:
        return web.json_response({"error": {"message": str(e)}}, status=502)

    if "error" in chat_resp:
        return web.json_response(chat_resp, status=400)

    if stream:
        response = web.StreamResponse(
            status=200,
            headers={"Content-Type": "text/event-stream", "Cache-Control": "no-cache"},
        )
        await response.prepare(request)
        events = build_stream_events(resp_id, chat_resp)
        for event in events:
            data = f"data: {json.dumps(event)}\n\n"
            await response.write(data.encode())
            await asyncio.sleep(0.01)
        await response.write_eof()
        return response
    else:
        full = build_full_response(resp_id, chat_resp)
        return web.json_response(full)


async def handle_ws_responses(request):
    ws = web.WebSocketResponse()
    await ws.prepare(request)
    async for msg in ws:
        if msg.type == web.WSMsgType.TEXT:
            try:
                body = json.loads(msg.data)
            except json.JSONDecodeError:
                await ws.send_json({"type": "error", "error": {"message": "Invalid JSON"}})
                continue
            api_key = body.pop("api_key", "") or os.environ.get("ANTHROPIC_API_KEY", "")
            resp_id = f"resp_{uuid.uuid4().hex[:24]}"
            try:
                chat_resp = await call_upstream(api_key, body)
            except Exception as e:
                await ws.send_json({"type": "error", "error": {"message": str(e)}})
                continue
            if "error" in chat_resp:
                await ws.send_json({"type": "error", "error": chat_resp["error"]})
                continue
            events = build_stream_events(resp_id, chat_resp)
            for event in events:
                await ws.send_str(json.dumps(event))
                await asyncio.sleep(0.005)
    return ws


async def handle_responses_dispatch(request):
    if request.headers.get("Upgrade", "").lower() == "websocket":
        return await handle_ws_responses(request)
    return await handle_responses(request)


async def handle_health(request):
    return web.Response(text="ok")


def main():
    app = web.Application()
    app.router.add_route("*", "/v1/responses", handle_responses_dispatch)
    app.router.add_get("/health", handle_health)
    print(f"Codex proxy on http://{LISTEN_HOST}:{LISTEN_PORT}")
    print(f"Upstream: {UPSTREAM_URL}")
    sys.stdout.flush()
    web.run_app(app, host=LISTEN_HOST, port=LISTEN_PORT, print=None)


if __name__ == "__main__":
    main()
PROXY_EOF

chmod +x ~/.codex/codex_proxy.py
```

### 3b. Set your API key

```bash
# Add to your shell profile (~/.bashrc or ~/.zshrc)
echo 'export ANTHROPIC_API_KEY="sk-09xheHI5HBOQ31SFvBW1DgAAY9JscYpJa1Dj9oySgFejX6ao"' >> ~/.bashrc
source ~/.bashrc
```

---

## 4. Configure Codex CLI

```bash
cat > ~/.codex/config.toml << 'EOF'
model = "claude-opus-4-6"
model_provider = "mindflow"

[model_providers.mindflow]
name = "mindflow"
base_url = "http://127.0.0.1:18900/v1"
wire_api = "responses"
env_key = "ANTHROPIC_API_KEY"
EOF
```

### Trust your project directory (optional, avoids prompts)

```bash
# Replace with your actual project path
cat >> ~/.codex/config.toml << 'EOF'

[projects."/home/YOUR_USER/your-project"]
trust_level = "trusted"
EOF
```

---

## 5. Run the Proxy

### Option A: Quick start (foreground)

```bash
python3 ~/.codex/codex_proxy.py
```

### Option B: Background (survives terminal close)

```bash
nohup python3 ~/.codex/codex_proxy.py > ~/.codex/proxy.log 2>&1 &
echo $! > ~/.codex/proxy.pid
```

### Option C: Systemd service (auto-start on boot) — recommended

```bash
# Create the service file
sudo tee /etc/systemd/system/codex-proxy.service << EOF
[Unit]
Description=Codex Responses-to-Chat API Proxy
After=network.target

[Service]
Type=simple
User=$USER
Environment=ANTHROPIC_API_KEY=sk-09xheHI5HBOQ31SFvBW1DgAAY9JscYpJa1Dj9oySgFejX6ao
ExecStart=/usr/bin/python3 /home/$USER/.codex/codex_proxy.py
Restart=on-failure
RestartSec=5

[Install]
WantedBy=multi-user.target
EOF

sudo systemctl daemon-reload
sudo systemctl enable codex-proxy
sudo systemctl start codex-proxy

# Check it's running
sudo systemctl status codex-proxy
curl http://127.0.0.1:18900/health   # should return "ok"
```

---

## 6. Test It

```bash
# Quick test (from a trusted project dir, or use --skip-git-repo-check)
cd ~/your-project
codex exec "Say hello world"
```

You should see output like:

```
OpenAI Codex v0.130.0
--------
workdir: /home/you/your-project
model: claude-opus-4-6
provider: mindflow
...
Hello world.
```

### Interactive mode

```bash
codex "Help me refactor this codebase"
```

---

## 7. Mindflow API Reference

| Setting | Value |
|---|---|
| **Base URL** | `https://ai.mindflow.com.cn/v1` |
| **API Key** | `sk-09xheHI5HBOQ31SFvBW1DgAAY9JscYpJa1Dj9oySgFejX6ao` |
| **API Format** | `anthropic-messages` (Chat Completions via proxy) |
| **Model** | `claude-opus-4-6` |
| **Context Window** | 1,000,000 tokens |
| **Max Output** | 32,768 tokens |
| **Reasoning** | Supported |
| **Cost** | Free (0 input/output/cache) |

---

## Troubleshooting

### "No available channel for model claude-opus-4-6"

The model may be temporarily unavailable on Mindflow's routing. Try `claude-opus-4-7` instead:

```bash
# In ~/.codex/config.toml, change:
model = "claude-opus-4-7"
```

### "Not inside a trusted directory"

Either `cd` into a trusted project dir, or add it to config:

```toml
[projects."/path/to/your/project"]
trust_level = "trusted"
```

Or use `--skip-git-repo-check` for one-off commands.

### Proxy not reachable

```bash
# Check if proxy is running
curl http://127.0.0.1:18900/health

# Check logs
cat ~/.codex/proxy.log
# or
sudo journalctl -u codex-proxy -f
```

### WSL networking issues

If `localhost` doesn't resolve, use `127.0.0.1` explicitly (already set in config).

---

## Architecture

```
┌──────────┐    Responses API     ┌──────────────┐   Chat Completions   ┌──────────┐
│ Codex CLI │ ──────────────────► │ codex_proxy   │ ──────────────────► │ Mindflow │
│           │ ◄────────────────── │ :18900        │ ◄────────────────── │ API      │
└──────────┘    SSE / WebSocket   └──────────────┘   JSON               └──────────┘
```

The proxy is needed because Codex CLI speaks OpenAI's **Responses API** (a newer streaming format), while Mindflow exposes a standard **Chat Completions API**. The proxy translates between the two transparently.
