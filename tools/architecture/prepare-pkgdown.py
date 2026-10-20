#!/usr/bin/env python3
"""Stage canonical developer Markdown as pkgdown website-only articles.

No package code is loaded or executed. Generated articles belong to the
existing vignettes/articles convention and are excluded from R distributions.
"""
from pathlib import Path
from urllib.parse import urlsplit, unquote, quote
import hashlib
import json
import posixpath
import re
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[2]
ARTICLES = ROOT / 'vignettes/articles'
ASSETS = ROOT / 'pkgdown/assets/contributor'
BASE = 'abdb5c713b089c3b10dea93bd94955f0f20c0988'
SPECIAL = {
    'tools/architecture/contributor-walkthrough.md': 'contributor-walkthrough',
    'tools/architecture/contributor-change-map.md': 'contributor-change-map',
    'tools/architecture/contributor-site-notes.md': 'contributor-site-notes',
    'tools/components/schema3/ARCHITECTURE.md': 'contributor-schema-architecture',
    'tools/components/schema3/OBJECT_FIELDS.md': 'contributor-object-fields',
    'tools/components/schema3/FIELD_DECISIONS.md': 'contributor-legacy-fields',
}
ROOT_DOCS = {'ARCHITECTURE.md', 'CONTRIBUTING.md', 'README.md',
             'CODE_OF_CONDUCT.md', 'LICENSE.md', 'SUPPORT.md',
             'ROADMAP.md', 'GOVERNANCE.md'}

def slug(path):
    rel = path.relative_to(ROOT).as_posix()
    return SPECIAL.get(rel, 'contributor-' + re.sub('[^a-z0-9]+', '-',
                       str(path.relative_to(ROOT).with_suffix('')).lower()).strip('-'))

queue = [ROOT / x for x in ['ARCHITECTURE.md', 'CONTRIBUTING.md', *SPECIAL]]
seen = set()
sources = {}
missing = []
previous_map = {}
prior_manifest = ROOT / 'tools/architecture/pkgdown-pages.json'
if prior_manifest.exists():
    previous_map = {Path(p['url']).stem: ROOT / p['source']
                    for p in json.loads(prior_manifest.read_text())['pages']}

def locate(origin, url):
    s = urlsplit(url)
    if s.scheme or s.netloc or not s.path:
        return None, s
    path = (origin.parent / unquote(s.path)).resolve()
    if not path.is_relative_to(ROOT):
        return None, s
    if path.is_dir() and (path / 'README.md').exists():
        path = path / 'README.md'
    if not path.exists() and path.suffix == '.html':
        # Links already normalized for the website resolve back to source.
        name = path.stem
        if path.parent == ROOT / 'articles':
            if name in previous_map:
                return previous_map[name], s
            for rel, stem in SPECIAL.items():
                if stem == name:
                    return ROOT / rel, s
            for old in seen | set(queue):
                if slug(old) == name:
                    return old, s
        for ext in ('.md', '.qmd'):
            if path.with_suffix(ext).exists():
                path = path.with_suffix(ext)
                break
    return path, s

while queue:
    path = queue.pop(0).resolve()
    if path in seen:
        continue
    seen.add(path)
    sources[path] = path.read_text()
    for url in re.findall(r'\]\(([^\s)]+)', sources[path]):
        p, s = locate(path, url)
        if p and p.is_file() and p.suffix.lower() in ('.md', '.qmd'):
            if p not in seen:
                queue.append(p)

def destination(path):
    rel = path.relative_to(ROOT).as_posix()
    if rel in ROOT_DOCS:
        return path.with_suffix('.html').name if path.name != 'README.md' else 'index.html'
    return 'articles/' + slug(path) + '.html'

def rewrite(text, origin, current):
    def link(match):
        prefix, url, suffix = match.groups()
        p, s = locate(origin, url)
        if p is None:
            return match.group(0)
        if p in sources:
            dst = destination(p)
        elif p.is_file() and p.suffix.lower() in ('.svg', '.png', '.jpg', '.jpeg', '.gif', '.webp'):
            name = re.sub('[^a-zA-Z0-9._-]+', '-', p.relative_to(ROOT).as_posix())
            shutil.copyfile(p, ASSETS / name)
            dst = 'contributor/' + name
        elif p.is_file() or p.is_dir():
            kind = 'tree' if p.is_dir() else 'blob'
            new = 'https://github.com/Katsevich-Lab/sceptre/' + kind + '/' + BASE + '/' + quote(p.relative_to(ROOT).as_posix())
            if s.fragment:
                new += '#' + s.fragment
            return prefix + new + suffix
        elif s.path.startswith(('articles/', 'contributor/')):
            # Explicit site destinations added to root Markdown.
            dst = s.path
        else:
            missing.append({'source': str(origin.relative_to(ROOT)), 'link': url})
            new = 'https://github.com/Katsevich-Lab/sceptre/blob/' + BASE + '/' + quote(p.relative_to(ROOT).as_posix())
            return prefix + new + ('#' + s.fragment if s.fragment else '') + suffix
        new = posixpath.relpath(dst, posixpath.dirname(current) or '.')
        if s.fragment:
            new += '#' + s.fragment
        return prefix + new + suffix
    return re.sub(r'(\]\()([^\s)]+)(\))', link, text)

ARTICLES.mkdir(parents=True, exist_ok=True)
ASSETS.mkdir(parents=True, exist_ok=True)
manifest = []
for path, text in sources.items():
    rel = path.relative_to(ROOT).as_posix()
    if rel in ROOT_DOCS:
        # Reuse existing top-level pkgdown pages. No second contribution portal.
        if rel in ('ARCHITECTURE.md', 'CONTRIBUTING.md'):
            path.write_text(rewrite(text, path, destination(path)))
        continue
    original = text
    if path.suffix == '.qmd':
        text = re.sub(r'\A---\n.*?\n---\n', '', text, flags=re.S)
        text = re.sub(r'```\{r[^}]*\}', '```r', text)
        text = re.sub(r'(?m)^#\|.*\n', '', text)
        text = re.sub(r'\{width=[^}]+\}', '', text)
    heading = re.search(r'^# (.+)$', text, re.M)
    title = heading.group(1) if heading else path.stem.replace('-', ' ').capitalize()
    if heading:
        text = text[:heading.start()] + text[heading.end():]
    text = rewrite(text.strip(), path, destination(path))
    historical = path.suffix == '.qmd' or any(x in rel.lower() for x in ('qualification', 'validation', 'review', 'disposition'))
    notice = ('> Supporting record: claims apply to the revisions and environments named below; this page does not qualify a later candidate.\n\n' if historical else '')
    if path.suffix == '.qmd':
        notice += '> Historical Quarto source rendered as static narrative and source blocks. Its R chunks are not executed.\n\n'
    nav = ('[Contributing](../CONTRIBUTING.html) · [Object model](../ARCHITECTURE.html#object-and-workflow-map) · [Discovery walkthrough](contributor-walkthrough.html) · [Object fields](contributor-object-fields.html)\n\n')
    header = ('---\ntitle: ' + json.dumps(title) + '\noutput: rmarkdown::html_document\n---\n\n'
              '<!-- Generated by tools/architecture/prepare-pkgdown.py from ' + rel + '; edit the canonical source. -->\n\n')
    out = ARTICLES / (slug(path) + '.Rmd')
    out.write_text(header + nav + notice + text + '\n')
    assert not re.search(r'```\{[rR]', out.read_text()), 'Executable chunk in generated article'
    manifest.append({'source': rel, 'source_sha256': hashlib.sha256(original.encode()).hexdigest(),
                     'article': str(out.relative_to(ROOT)), 'url': destination(path)})

(ROOT / 'tools/architecture/pkgdown-pages.json').write_text(json.dumps({'source_base': BASE, 'pages': manifest, 'missing_original_targets': missing}, indent=2) + '\n')
print(json.dumps({'generated_articles': len(manifest), 'missing_original_targets': missing}, indent=2))
