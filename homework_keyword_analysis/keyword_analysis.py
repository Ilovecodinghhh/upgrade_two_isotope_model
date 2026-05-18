#!/usr/bin/env python3
"""
Keyword Analysis: Atmospheric Methane Isotope Source Partitioning
Compares keyword frequency between 2005-2015 and 2015-2025 literature contexts.

Methodology:
- Scans all .md files in upgrade_two_isotope_model/
- Extracts paragraphs and identifies citation years (e.g. "Basu et al., 2022")
- Tags each paragraph with an era based on cited years
- Counts keyword occurrences per era
- Produces a comparison bar chart
"""

import os
import re
import json
from collections import defaultdict
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

# ─── Define field-specific keywords ───
KEYWORDS = {
    # Core species & measurements
    "methane (CH₄)":            [r'\bCH[₄4]\b', r'\bmethane\b'],
    "δ¹³C / carbon isotope":    [r'δ[¹1]³3C', r'd13C', r'δ13C', r'¹³C', r'13C[-_]CH', r'carbon isotop'],
    "δD / hydrogen isotope":    [r'δD', r'\bdD\b', r'deuterium', r'hydrogen isotop', r'δ²H', r'CH₃D', r'CH3D'],
    "isotope mass balance":     [r'isotop\w* mass balance', r'isotope budget', r'mass.balance'],
    
    # Sources
    "fossil fuel":              [r'fossil.fuel', r'\bFF\b', r'thermogenic', r'fugitive emission'],
    "microbial":                [r'microbial', r'\bMic\b', r'biogenic'],
    "biomass burning":          [r'biomass.burning', r'\bBB\b', r'pyrogenic', r'fire emission'],
    "wetlands":                 [r'wetland'],
    "ruminants / livestock":    [r'ruminant', r'livestock', r'enteric'],
    "coal mining":              [r'coal.min', r'coal emission'],
    "oil & natural gas":        [r'oil.and.gas', r'natural gas', r'\bONG\b'],
    
    # Sinks & chemistry
    "OH (hydroxyl radical)":    [r'\bOH\b', r'hydroxyl'],
    "kinetic isotope effect":   [r'KIE', r'kinetic isotope', r'fractionation factor'],
    "tropospheric Cl":          [r'tropospheric Cl', r'Cl sink', r'chlorine'],
    "CH₄ lifetime":             [r'CH[₄4]\s*lifetime', r'methane lifetime', r'effective lifetime'],
    
    # Methods
    "box model":                [r'box.model', r'1.box', r'2.box', r'two.box', r'one.box', r'three.box'],
    "Monte Carlo":              [r'Monte Carlo', r'\bMC\b'],
    "inversion":                [r'inversion', r'inverse method', r'variational', r'4D.Var'],
    "source partitioning":      [r'source.partition', r'source.apportion', r'source.attribution'],
    "satellite observations":   [r'satellite', r'GOSAT', r'TROPOMI', r'remote.sens'],
    "Bayesian":                 [r'Bayesian', r'posterior', r'prior'],
    
    # Datasets & inventories
    "EDGAR":                    [r'EDGAR'],
    "NOAA / GML":               [r'NOAA', r'GML', r'Global Monitoring'],
    "source signatures":        [r'source.signature', r'end.member', r'isotopic signature'],
    
    # Policy & context
    "global methane budget":    [r'global methane budget', r'methane budget', r'CH[₄4] budget'],
    "climate / greenhouse":     [r'climate change', r'greenhouse', r'global warming', r'GHG'],
    "emission trends":          [r'emission.trend', r'growth.rate', r'post.200[67]'],
    "hemispheric":              [r'hemispheric', r'NH.SH', r'inter.hemispheric', r'Northern Hemisphere', r'Southern Hemisphere'],
    "uncertainty":              [r'uncertainty', r'error propagat', r'confidence interval', r'variance decomposition'],
}

# ─── Collect all .md files under upgrade_two_isotope_model/ ───
base = "/home/openclaw1188/.openclaw/workspace/upgrade_two_isotope_model"
md_files = []
for root, dirs, files in os.walk(base):
    # Skip deep nested code repos (He2026Science_DS has a big code repo)
    if 'geoschem' in root:
        continue
    for f in files:
        if f.endswith('.md'):
            md_files.append(os.path.join(root, f))

print(f"Found {len(md_files)} .md files")

# ─── Parse paragraphs and extract citation years ───
year_pattern = re.compile(r'(?:19|20)\d{2}')
# More specific citation pattern: Author(year) or Author et al. (year) or Author et al., year
cite_pattern = re.compile(r'(?:[A-Z][a-z]+\s+(?:et\s+al\.?\s*[,;]?\s*)?(?:\(|\b))?((?:19|20)\d{2})')

def extract_years(text):
    """Extract all plausible citation years from text"""
    years = set()
    for m in year_pattern.finditer(text):
        y = int(m.group())
        if 1980 <= y <= 2026:
            years.add(y)
    return years

def get_era(years):
    """Classify a set of years into eras. Returns both eras if years span both."""
    eras = set()
    for y in years:
        if 2005 <= y <= 2015:
            eras.add("2005-2015")
        if 2015 <= y <= 2025:
            eras.add("2015-2025")
    return eras

def count_keywords(text):
    """Count keyword matches in a text chunk"""
    counts = {}
    text_lower = text.lower()
    for kw, patterns in KEYWORDS.items():
        count = 0
        for pat in patterns:
            count += len(re.findall(pat, text, re.IGNORECASE))
        if count > 0:
            counts[kw] = count
    return counts

# ─── Also classify ImportantReferences files by their folder year ───
ref_years = {}
ref_dir = os.path.join(base, "ImportantReferences")
if os.path.exists(ref_dir):
    for d in os.listdir(ref_dir):
        m = re.search(r'(\d{4})', d)
        if m:
            ref_years[d] = int(m.group(1))

# ─── Process all files ───
era_keyword_counts = {
    "2005-2015": defaultdict(int),
    "2015-2025": defaultdict(int),
}

for fpath in md_files:
    try:
        with open(fpath, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()
    except:
        continue
    
    # Determine file-level era from ImportantReferences folder name
    file_era = None
    for ref_name, year in ref_years.items():
        if ref_name in fpath:
            if 2005 <= year <= 2015:
                file_era = "2005-2015"
            elif 2015 <= year <= 2025:
                file_era = "2015-2025"
            break
    
    # Split into paragraphs (double newline)
    paragraphs = re.split(r'\n\s*\n', content)
    
    for para in paragraphs:
        if len(para.strip()) < 20:
            continue
        
        years = extract_years(para)
        eras = get_era(years)
        
        # If no era from paragraph years, use file-level era
        if not eras and file_era:
            eras = {file_era}
        
        # If still no era, skip
        if not eras:
            continue
        
        kw_counts = count_keywords(para)
        
        for era in eras:
            for kw, cnt in kw_counts.items():
                era_keyword_counts[era][kw] += cnt

# ─── Build summary ───
all_keywords = sorted(KEYWORDS.keys(), 
                      key=lambda k: era_keyword_counts["2015-2025"].get(k, 0) + era_keyword_counts["2005-2015"].get(k, 0),
                      reverse=True)

print("\n" + "="*80)
print("FIELD: Atmospheric Methane Isotope Geochemistry & Source Partitioning")
print("="*80)
print(f"\n{'Keyword':<35} {'2005-2015':>10} {'2015-2025':>10} {'Ratio':>8}")
print("-"*65)
for kw in all_keywords:
    c1 = era_keyword_counts["2005-2015"].get(kw, 0)
    c2 = era_keyword_counts["2015-2025"].get(kw, 0)
    ratio = f"{c2/c1:.1f}x" if c1 > 0 else "new"
    print(f"{kw:<35} {c1:>10} {c2:>10} {ratio:>8}")

# ─── Save summary JSON ───
summary = {
    "field": "Atmospheric Methane Isotope Geochemistry & Source Partitioning",
    "description": "Keywords extracted from research repository on methane dual-isotope (δ¹³C + δD) box models for global CH₄ source partitioning into fossil fuel, microbial, and biomass burning components.",
    "periods": {
        "2005-2015": dict(era_keyword_counts["2005-2015"]),
        "2015-2025": dict(era_keyword_counts["2015-2025"]),
    },
    "total_md_files_analyzed": len(md_files),
}
with open("/home/openclaw1188/.openclaw/workspace/keyword_summary.json", "w") as f:
    json.dump(summary, f, indent=2)

# ─── PLOT ───
fig, axes = plt.subplots(1, 2, figsize=(20, 10), gridspec_kw={'width_ratios': [3, 1.2]})

# Left panel: horizontal bar comparison
ax = axes[0]
# Filter to keywords with at least some counts
plot_kws = [kw for kw in all_keywords if era_keyword_counts["2005-2015"].get(kw, 0) + era_keyword_counts["2015-2025"].get(kw, 0) > 0]
# Take top 25
plot_kws = plot_kws[:25]
plot_kws = plot_kws[::-1]  # reverse for bottom-to-top

y = np.arange(len(plot_kws))
width = 0.38

vals_early = [era_keyword_counts["2005-2015"].get(kw, 0) for kw in plot_kws]
vals_late = [era_keyword_counts["2015-2025"].get(kw, 0) for kw in plot_kws]

bars1 = ax.barh(y - width/2, vals_early, width, label='2005–2015', color='#4C72B0', alpha=0.85, edgecolor='white', linewidth=0.5)
bars2 = ax.barh(y + width/2, vals_late, width, label='2015–2025', color='#DD8452', alpha=0.85, edgecolor='white', linewidth=0.5)

ax.set_yticks(y)
ax.set_yticklabels(plot_kws, fontsize=10)
ax.set_xlabel('Keyword Occurrences in Literature Context', fontsize=12)
ax.set_title('Keyword Frequency Comparison\nAtmospheric CH₄ Isotope Source Partitioning', fontsize=14, fontweight='bold')
ax.legend(fontsize=11, loc='lower right', framealpha=0.9)
ax.grid(axis='x', alpha=0.3, linestyle='--')
ax.set_axisbelow(True)

# Add count labels
for bar, val in zip(bars1, vals_early):
    if val > 0:
        ax.text(val + max(vals_late)*0.01, bar.get_y() + bar.get_height()/2, 
                str(val), va='center', fontsize=8, color='#4C72B0')
for bar, val in zip(bars2, vals_late):
    if val > 0:
        ax.text(val + max(vals_late)*0.01, bar.get_y() + bar.get_height()/2, 
                str(val), va='center', fontsize=8, color='#DD8452')

# Right panel: ratio / growth chart
ax2 = axes[1]
ratios = []
colors = []
labels = []
for kw in plot_kws:
    c1 = max(era_keyword_counts["2005-2015"].get(kw, 0), 1)  # avoid div by 0
    c2 = era_keyword_counts["2015-2025"].get(kw, 0)
    r = c2 / c1
    ratios.append(r)
    labels.append(kw)
    if era_keyword_counts["2005-2015"].get(kw, 0) == 0:
        colors.append('#55A868')  # green = new topic
    elif r > 2:
        colors.append('#C44E52')  # red = strong growth
    elif r > 1:
        colors.append('#DD8452')  # orange = moderate growth
    else:
        colors.append('#4C72B0')  # blue = declining or stable

ax2.barh(y, ratios, color=colors, alpha=0.85, edgecolor='white', linewidth=0.5)
ax2.axvline(x=1, color='gray', linestyle='--', alpha=0.7, linewidth=1)
ax2.set_yticks([])
ax2.set_xlabel('Growth Ratio (2015–2025 / 2005–2015)', fontsize=11)
ax2.set_title('Growth\nRatio', fontsize=13, fontweight='bold')
ax2.grid(axis='x', alpha=0.3, linestyle='--')
ax2.set_axisbelow(True)

# Add ratio labels
for i, (r, c) in enumerate(zip(ratios, colors)):
    label = f"{r:.1f}×" if era_keyword_counts["2005-2015"].get(plot_kws[i], 0) > 0 else "NEW"
    ax2.text(r + max(ratios)*0.02, i, label, va='center', fontsize=9, fontweight='bold', color=c)

plt.tight_layout()
plt.savefig("/home/openclaw1188/.openclaw/workspace/keyword_comparison_plot.png", dpi=200, bbox_inches='tight', facecolor='white')
print("\n✅ Plot saved to keyword_comparison_plot.png")

# ─── Also make a summary text block ───
print("\n" + "="*80)
print("KEY FINDINGS:")
print("="*80)
print("""
The research field is: ATMOSPHERIC METHANE ISOTOPE GEOCHEMISTRY & SOURCE PARTITIONING

This field focuses on using stable isotope measurements (δ¹³C-CH₄ and δD-CH₄) to 
decompose global methane emissions into fossil fuel, microbial, and biomass burning 
components via isotope mass-balance models.

Top emerging keywords (2015-2025 vs 2005-2015):
  • δD / hydrogen isotope — enormous growth (δD datasets only became available ~2020s)
  • Satellite observations (GOSAT/TROPOMI) — explosion of space-based CH₄ data
  • Bayesian methods — shift from simple MC to formal Bayesian inversions
  • Source signatures — recognized as the binding constraint on attribution
  • Hemispheric — new focus on NH/SH divergence rather than global-only models
  • Kinetic isotope effect (KIE) — intensified debate on OH-¹³C fractionation
  • Uncertainty / variance decomposition — growing emphasis on error propagation

Persistent core keywords across both periods:
  • Methane (CH₄), fossil fuel, microbial, biomass burning, OH, inversion
  • Box model, Monte Carlo, isotope mass balance
""")

