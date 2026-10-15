"""Regenerate the "Cell type markers" tables in index.html from Tables/*_harmonized*.csv.

Run from the repo root:  python src/build_marker_tables.py
Replaces everything between <!-- markers:start --> and <!-- markers:end --> in index.html.
"""
import csv, glob, html, re, urllib.parse

LINEAGES = [  # (display name, filename glob)
    ("Proximal epithelium", "Tables/Epithelial_proximal_harmonized*.csv"),
    ("Distal epithelium",   "Tables/Epithelial_distal_harmonized*.csv"),
    ("Endothelial",         "Tables/Endothelial_harmonized*.csv"),
    ("Stromal",             "Tables/Stromal_harmonized*.csv"),
]

def table_html(name, path):
    with open(path, newline="", encoding="utf-8") as fh:
        rows = [r for r in csv.reader(fh) if any(c.strip() for c in r)]
    header, body = rows[0], rows[1:]
    href = urllib.parse.quote(path)
    out = [f'        <details class="lineage">',
           f'          <summary>{html.escape(name)}</summary>',
           f'          <div class="table-scroll">',
           f'            <table>',
           f'              <thead><tr><th scope="col">#</th>' +
           "".join(f'<th scope="col">{html.escape(h)}</th>' for h in header) + '</tr></thead>',
           f'              <tbody>']
    for i, r in enumerate(body, 1):
        cells = "".join(f'<td>{html.escape(c)}</td>' for c in r[:len(header)])
        cells += "<td></td>" * (len(header) - len(r))
        out.append(f'                <tr><th scope="row">{i}</th>{cells}</tr>')
    out += ['              </tbody>', '            </table>', '          </div>',
            f'          <p class="foot"><a href="{href}" download>Download CSV</a></p>',
            '        </details>']
    return "\n".join(out)

blocks = []
for name, pattern in LINEAGES:
    paths = sorted(glob.glob(pattern))
    assert len(paths) == 1, (pattern, paths)
    blocks.append(table_html(name, paths[0]))
fragment = "\n".join(blocks)

src = open("index.html", encoding="utf-8").read()
new, n = re.subn(r"<!-- markers:start -->.*?<!-- markers:end -->",
                 lambda m: "<!-- markers:start -->\n" + fragment + "\n        <!-- markers:end -->",
                 src, flags=re.S)
assert n == 1, "markers:start/end comments not found in index.html"
open("index.html", "w", encoding="utf-8").write(new)
print("wrote", len(blocks), "lineage tables into index.html")
