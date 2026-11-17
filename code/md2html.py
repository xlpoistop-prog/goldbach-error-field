import re, html, sys

def inline(s):
    s = html.escape(s)
    # protect code spans first: emphasis must never look inside them
    spans = []

    def _keep(m):
        spans.append(m.group(1))
        return "\x00%d\x00" % (len(spans) - 1)

    s = re.sub(r'`([^`]+)`', _keep, s)
    s = re.sub(r'\*\*([^*]+)\*\*', r'<strong>\1</strong>', s)
    s = re.sub(r'(?<!\*)\*([^*]+)\*(?!\*)', r'<em>\1</em>', s)
    s = re.sub(r'\x00(\d+)\x00', lambda m: '<code>' + spans[int(m.group(1))] + '</code>', s)
    return s

def convert(md):
    lines = md.split('\n')
    out = []; i = 0
    while i < len(lines):
        ln = lines[i]
        if ln.startswith('```'):
            i += 1; buf = []
            while i < len(lines) and not lines[i].startswith('```'):
                buf.append(html.escape(lines[i])); i += 1
            out.append('<pre>' + '\n'.join(buf) + '</pre>'); i += 1
            continue
        if ln.startswith('|') and i+1 < len(lines) and set(lines[i+1].replace('|','').replace('-','').replace(':','').strip()) <= set():
            hdr = [c.strip() for c in ln.strip().strip('|').split('|')]
            i += 2; rows = []
            while i < len(lines) and lines[i].startswith('|'):
                rows.append([c.strip() for c in lines[i].strip().strip('|').split('|')]); i += 1
            t = '<table><thead><tr>' + ''.join('<th>'+inline(h)+'</th>' for h in hdr) + '</tr></thead><tbody>'
            for r in rows:
                t += '<tr>' + ''.join('<td>'+inline(c)+'</td>' for c in r) + '</tr>'
            out.append(t + '</tbody></table>')
            continue
        mimg = re.match(r'^!\[(.*?)\]\((.*?)\)$', ln.strip())
        if mimg:
            out.append('<figure><img src="%s" alt="%s"/><figcaption>%s</figcaption></figure>'
                       % (html.escape(mimg.group(2)), html.escape(mimg.group(1)), inline(mimg.group(1))))
            i += 1; continue
        m = re.match(r'^(#{1,4})\s+(.*)$', ln)
        if m:
            lv = len(m.group(1)); out.append(f'<h{lv}>{inline(m.group(2))}</h{lv}>'); i += 1; continue
        if ln.strip() == '---':
            out.append('<hr/>'); i += 1; continue
        if ln.strip() == '':
            i += 1; continue
        if re.match(r'^\s*[-*]\s+', ln) or re.match(r'^\s*\d+\.\s+', ln):
            items = []
            bullet = lambda t: bool(re.match(r'^\s*[-*]\s+', t) or re.match(r'^\s*\d+\.\s+', t))
            while i < len(lines) and bullet(lines[i]):
                it = re.sub(r'^\s*(?:[-*]|\d+\.)\s+', '', lines[i]); i += 1
                while (i < len(lines) and lines[i].strip() and not bullet(lines[i])
                       and not re.match(r'^(#{1,4}\s|\||```|---)', lines[i])):
                    it += ' ' + lines[i].strip(); i += 1
                items.append(it)
            out.append('<ul>' + ''.join('<li>'+inline(x)+'</li>' for x in items) + '</ul>'); continue
        para = [ln]; i += 1
        while i < len(lines) and lines[i].strip() and not re.match(r'^(#{1,4}\s|\||```|---)', lines[i]):
            para.append(lines[i]); i += 1
        out.append('<p>' + inline(' '.join(para)) + '</p>')
    return '\n'.join(out)

CSS = """

body{max-width:860px;margin:2em auto;padding:0 1.2em;font-family:Georgia,"Songti SC",serif;
line-height:1.75;color:#1a1a1a;background:#fbfbf9}
h1{font-size:1.7em;border-bottom:2px solid #333;padding-bottom:.3em;margin-top:1.6em}
h2{font-size:1.35em;margin-top:1.8em;color:#222;border-bottom:1px solid #ddd;padding-bottom:.2em}
h3{font-size:1.12em;margin-top:1.4em;color:#333}
h4{font-size:1.02em;margin-top:1.2em;color:#444}
pre{background:#f4f4f0;border-left:3px solid #bbb;padding:.7em 1em;overflow-x:auto;
font-family:"SF Mono",Menlo,monospace;font-size:.86em;line-height:1.5}
code{font-family:"SF Mono",Menlo,monospace;font-size:.9em;background:#f2f2ee;padding:.08em .25em;border-radius:3px}
pre code{background:none;padding:0}
table{border-collapse:collapse;margin:1em 0;font-size:.92em;width:100%}
th,td{border:1px solid #ccc;padding:.35em .6em;text-align:left}
th{background:#eee;font-weight:600}
tbody tr:nth-child(even){background:#f7f7f4}
hr{border:none;border-top:1px solid #ddd;margin:2em 0}
blockquote{border-left:3px solid #bbb;margin:1em 0;padding:.2em 1em;color:#555;background:#f6f6f2}
ul{margin:.6em 0 .6em 1.2em}
p{margin:.8em 0}
figure{margin:1.4em 0;text-align:center}
figure img{max-width:100%;border:1px solid #e2e2dd;background:#fff}
figcaption{font-size:.88em;color:#555;margin-top:.4em}
@media print{
  @page{size:A4;margin:20mm 18mm 18mm 18mm}
  body{max-width:none;margin:0;padding:0;background:#fff;color:#000;
       font-family:"Times New Roman","Songti SC",Georgia,serif;font-size:10.5pt;line-height:1.55}
  h1{font-size:16pt;border-bottom:1.2pt solid #000;margin:0 0 .6em 0;page-break-after:avoid}
  h2{font-size:12.5pt;border-bottom:.6pt solid #999;page-break-after:avoid;margin-top:1.4em}
  h3{font-size:11pt;page-break-after:avoid}
  h4{font-size:10.5pt;page-break-after:avoid}
  p{margin:.55em 0;text-align:justify;orphans:2;widows:2}
  ul{margin:.4em 0 .4em 1.1em}
  table{font-size:8.6pt;page-break-inside:avoid;margin:.8em 0}
  th,td{padding:.22em .4em;border:.5pt solid #666}
  th{background:#eee !important;-webkit-print-color-adjust:exact}
  tbody tr:nth-child(even){background:#f8f8f8 !important;-webkit-print-color-adjust:exact}
  pre,code{font-family:Consolas,"DejaVu Sans Mono","Cambria Math",monospace}
  pre{font-size:8.8pt;background:none !important;border:none;border-left:1.2pt solid #ccc;
      padding:.25em .0em .25em .9em;margin:.5em 0 .6em .4em;page-break-inside:avoid;
      white-space:pre-wrap;word-break:break-word}
  code{font-size:9.6pt;background:none !important;padding:0;border-radius:0}
  figure{page-break-inside:avoid;margin:1em 0}
  figure img{max-width:78%;border:.5pt solid #bbb}
  figcaption{font-size:8.6pt}
  hr{display:none}
  a{color:#000;text-decoration:none}
}
"""

md = open(sys.argv[1], encoding='utf-8').read()
title = sys.argv[3] if len(sys.argv) > 3 else '哥德巴赫误差场的乘性结构'
lang = sys.argv[4] if len(sys.argv) > 4 else 'zh-CN'
html_out = f"""<!DOCTYPE html>
<html lang="{lang}"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{html.escape(title)}</title><style>{CSS}</style></head><body>
{convert(md)}
</body></html>"""
open(sys.argv[2], 'w', encoding='utf-8').write(html_out)
print("written", sys.argv[2], len(html_out))
