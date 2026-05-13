# Codex CLI — Third-Party API Setup Guide

Connect OpenAI's Codex CLI to any OpenAI-compatible API provider (e.g., MindFlow, DeepSeek, Qwen, etc.).

---

## Method 1: Config File (Recommended)

### 1. Create/edit the config file

Path: `~/.codex/config.toml`

```toml
model_provider = "mindflow"          # Custom provider name
model = "gpt-5.4"                    # Model to use
model_reasoning_effort = "medium"    # low / medium / high

[model_providers.mindflow]
name = "mindflow"                    # Must match model_provider above
base_url = "https://api.mindflow.example/v1"  # Provider endpoint (include /v1)
env_key = "MINDFLOW_API_KEY"         # Name of the env var holding your key
```

### 2. Set your API key

```bash
# Linux / macOS
export MINDFLOW_API_KEY="sk-your-key-here"

# Windows CMD
set MINDFLOW_API_KEY=sk-your-key-here

# Windows PowerShell
$env:MINDFLOW_API_KEY="sk-your-key-here"
```

> **Tip:** Add the export line to `~/.bashrc` or `~/.zshrc` to make it persistent.

### 3. Verify & restart

```bash
echo $MINDFLOW_API_KEY   # Should print your key
```

Fully quit and reopen VS Code (or your terminal) for changes to take effect.

---

## Method 2: Open-Codex Fork (Advanced)

The community fork [open-codex](https://github.com/open-codex/open-codex) replaces the Responses API with the standard **Chat Completions API**, enabling compatibility with any OpenAI-style provider.

```bash
npm install -g open-codex
```

Config is the same as Method 1, but add `wire_api = "chat_completions"` inside your provider block:

```toml
[model_providers.mindflow]
name = "mindflow"
base_url = "https://api.mindflow.example/v1"
env_key = "MINDFLOW_API_KEY"
wire_api = "chat_completions"
```

---

## Troubleshooting

| Error | Likely Cause | Fix |
|-------|-------------|-----|
| **401** | Invalid API key | Check for typos, extra spaces, or newlines in key |
| **403** | No access to model / no balance | Verify model name and account balance |
| **503** | Wrong endpoint | Ensure `base_url` ends with `/v1` |

### Other tips

- **Free-tier limits:** Codex makes frequent tool calls — free quotas burn fast. Use sparingly or switch to pay-per-token plans.
- **Visual config tools:** Projects like CC-Switch offer a GUI for managing API configs across multiple AI tools.

---

## Quick Install Summary

```bash
# 1. Install Codex CLI
npm install -g @openai/codex

# 2. Create config
mkdir -p ~/.codex
cat > ~/.codex/config.toml << 'EOF'
model_provider = "mindflow"
model = "gpt-5.4"
model_reasoning_effort = "medium"

[model_providers.mindflow]
name = "mindflow"
base_url = "https://api.mindflow.example/v1"
env_key = "MINDFLOW_API_KEY"
EOF

# 3. Set API key
export MINDFLOW_API_KEY="sk-your-key-here"

# 4. Run
codex
```
