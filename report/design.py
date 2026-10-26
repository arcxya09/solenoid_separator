# 设计报告：复用 report_ex2.html 的样式和绘图脚本，正文用 design_body.html
import re, sys
src = open('report_ex2.html').read()
body = open(sys.argv[1] if len(sys.argv)>1 else 'design_body.html').read()
i = src.index('<div class="wrap">'); j = src.index('<script', i)
head, oldbody, scripts = src[:i], src[i:j], src[j:]
head = head.replace('<title>螺线管同位素分离器</title>', '<title>螺线管分离器设计报告</title>')
assert '设计报告' in head
lbl = 'Rb-87 减速口袋（初版 300 eV 已撤回，改 ≥ 1.5 keV），冷底面'
assert scripts.count(lbl) == 1
scripts = scripts.replace(lbl, 'Rb-87 产品口袋（落点 ≥ 1.5 keV），冷底面')
# 正文未用到、但脚本会访问的元素放进隐藏容器
ids = set(re.findall(r'id="([^"]+)"', oldbody)) - set(re.findall(r'id="([^"]+)"', body))
hidden = '<div hidden aria-hidden="true" style="display:none">' + ''.join(
    f'<table><tbody id="{k}"></tbody></table>' if k in ('tol-body', 'f-hv', 'h-cx') else
    f'<svg id="{k}"></svg>' if k in ('colsvg','schem') else f'<div id="{k}" class="plot"></div>'
    for k in sorted(ids) if not re.fullmatch(r's\d+', k)) + '</div>\n'
open(sys.argv[2] if len(sys.argv)>2 else 'report_design.html', 'w').write(head + body + hidden + scripts)
print('hidden ids:', sorted(ids))
