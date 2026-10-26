import re, sys
M = open('main.html').read().split('\n')
B = open('bl.html').read().split('\n')

def R(src, a, b):
    """lines a..b inclusive, 1-indexed"""
    return '\n'.join(src[a-1:b])

def inner(src, s, e):
    """section body without <section>, secno and h2 lines (s=<section> line, e=</section> line)"""
    return R(src, s+3, e-1)

def sub(t, old, new, count=1):
    if old not in t:
        sys.exit('MISSING: ' + old[:80])
    return t.replace(old, new, count)

def sec(no, title, h2, body, sid=None):
    i = f' id="{sid}"' if sid else ''
    return f'<section{i}>\n  <div class="secno">{no} · {title}</div>\n  <h2>{h2}</h2>\n{body}\n</section>\n'

# ---------------- head ----------------
head_end = next(i for i, l in enumerate(M) if l.startswith('<div class="wrap">'))
head = '\n'.join(M[1:head_end])
head = sub(head, '<title>螺线管同位素分离器</title>', '<title>螺线管同位素分离器</title>')
# beamline extra css block
bi = next(i for i, l in enumerate(B) if l.startswith('/* 本页补充 */'))
bj = next(i for i in range(bi, len(B)) if B[i].startswith('</style>'))
extra_css = '\n'.join(B[bi:bj]) + '\n.toc{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:2px 24px;font-size:14px;margin:0;padding-left:1.4em;max-width:76ch}\n@media (max-width:640px){.toc{grid-template-columns:1fr}}\n.upd{font-size:12px;font-family:var(--f-mono);color:var(--warn);margin-left:6px}\n/* 文字与图同宽 */\np,.lede,ul.tight,ol.steps,.note,.caution,figcaption,.toc{max-width:none}'
extra_css = extra_css.replace('.plot svg{', '.plot > svg{')
head = sub(head, '</style>', extra_css + '\n</style>')

# ---------------- header + summary (new) ----------------
header = '''<div class="wrap">
<header>
  <div class="eyebrow">项目报告 · 概念设计与数值验证 · Rb-85 / Rb-87</div>
  <h1>用超导螺线管取代大型二极磁铁的电磁同位素分离器</h1>
  <p class="lede">思路来自你的原始设想：离子束在螺线管中回旋，不同质量转过的相位略有不同。分析和全三维粒子追踪表明，这一效应确实可以用来分离同位素：采用<b>环形平行束</b>时，螺线管变成一个焦距与质量相关的轴对称透镜：Rb-87 汇聚到轴上，Rb-85 落在外围圆环上，用一个中心孔即可单级取出高丰度 Rb-87。本报告汇总了原理、离子源、束流线、收集器、空间电荷、真空和流强放大的全部结果，全篇以当前基准设计为准。</p>
  <div class="kpis">
    <div class="kpi"><b>≥ 99.97%</b><span>Rb-87 产品丰度（天然 27.8%，单级，未计散射）</span></div>
    <div class="kpi"><b>≈ 100%</b><span>Rb-87 收集效率（中心孔半径 5.9 mm）</span></div>
    <div class="kpi"><b>2.785 T × 1 m</b><span>超导螺线管，室温孔 ≥ 180 mm，Rb⁺ 10 keV</span></div>
    <div class="kpi"><b>0.3 / 3.5 mA</b><span>无中和流强上限：基准 / 放大设计</span></div>
  </div>
</header>
'''

summary = sec('00', '摘要', '要点', '''  <ul class="tight">
    <li><b>原理成立。</b>从无场区进入螺线管的离子，其横向运动在以 Larmor 频率旋转的坐标系中是简谐振荡，总转角 α = q∫B dz /(2p) ∝ m<sup>−1/2</sup>。对 Rb-85/87，Δα/α = Δm/2m ≈ 1.16%，基准设计中 Δα ≈ 0.12 rad（第 01 节）。</li>
    <li><b>环形平行束。</b>所有离子从同一半径平行入射，螺旋运动带来的非近轴相位对每个离子都相同，不形成像差，螺线管成为一个焦距与质量相关的透镜（第 02 节）。</li>
    <li><b>基准设计：</b>多孔钨（Re/Ir 表面）环形表面电离源，R = 30 mm，+10 kV；源出口 z = −0.6 m；超导螺线管 2.785 T、长 1 m；收集器在 z = 1.365 m，中心孔半径 5.9 mm。焦面上 Rb-87 最大半径 5.5 mm，Rb-85 最小半径 7.4 mm，Rb-87 收集效率约 100%，丰度 ≥ 99.97%（第 03、07 节）。</li>
    <li><b>离子源：</b>Pierce 电极加环形单透镜（−16 kV）给出径向/方位角散 1.3 / 2.5 mrad、发射度 0.83 mm·mrad。工程核算没有发现原理障碍：Re/Ir 表面 1400 K 时电离效率 ≥ 99.7%，电极倒圆角后最大表面场 3.1 kV/mm，加热约 75 W，0.3 mA 时耗铷 23 mg/天（第 04、05 节）。</li>
    <li><b>束流线很短、很简单。</b>源出口到线圈 0.6 m 漂移，不需要四极或匹配透镜。源区边缘场（16 mT）不需要屏蔽，发射面残余场 −10 到 +25 mT 都没有影响。两组静电导向器加可平移收集器可校正线圈 5 mm 偏心或 10 mrad 倾斜（第 06 节）。</li>
    <li><b>最紧的容差</b>是线圈对准（不加导向时偏心 ≤ 1 mm、倾斜 ≤ 1 mrad）、源倾斜（≤ 3 mrad）、B₀ ±0.18% 和源高压 ±0.2%。各项主要通过焦点位置耦合，共享收集器 −20 到 +10 mm 的纵向容差（第 08 节）。</li>
    <li><b>流强是最紧的约束。</b>无中和时基准设计上限约 0.3 mA；焦点对未中和流强的灵敏度约 0.9 m/mA，收集器固定时只能容许 ±0.05 mA 的波动。因此中和不仅要做到，还要<b>稳定</b>：30 mA 总束流时中和度需稳定在 ±0.17% 以内（第 09 节）。</li>
    <li><b>丰度的主要威胁是热铷原子，不是散射。</b>电荷交换产生的 Rb-85 快中性原子不会进入收集孔；大角散射在 10⁻⁵ Pa 时贡献约 5×10⁻⁶。收集区天然铷蒸气直接沉积才是主要风险，收集区和废束靶需冷到 ≤ 230–250 K（第 10 节）。</li>
    <li><b>流强放大：</b>无中和上限按 I ∝ E<sup>1.5</sup>R² 增长。放大几何（30 keV、R = 50 mm、1.7 m × 2.8 T）后约 3.5 mA，已模拟验证；保持几何只提场强（40 keV、5.6 T）约 2.4 mA。达到 10 / 100 mA 分别需要约 70% / 96% 的中和，低于 calutron 实际达到的约 99.7%。50 mA 时 Rb-87 产量约 1 g/天（第 11 节）。</li>
    <li><b>唯一需要实验确认的关键未知数</b>是 3 T 强磁场中电子中和的效果和稳定度，建议作为第一项实验（第 13 节）。</li>
  </ul>
  <h3>目录</h3>
  <ol class="toc" start="0">
    <li><a href="#s00">摘要</a></li><li><a href="#s01">原理</a></li><li><a href="#s02">方案概念与第一轮结果</a></li><li><a href="#s03">整机布局与基准参数</a></li>
    <li><a href="#s04">离子源：引出光学</a></li><li><a href="#s05">离子源：工程可行性</a></li><li><a href="#s06">束流线：边缘场、调谐与对准</a></li><li><a href="#s07">收集器与分离性能</a></li>
    <li><a href="#s08">容差汇总</a></li><li><a href="#s09">空间电荷与中和</a></li><li><a href="#s10">真空、电荷交换与铷蒸气</a></li><li><a href="#s11">流强放大与产量</a></li>
    <li><a href="#s13">与传统分离器对比、成本</a></li><li><a href="#s14">风险与下一步</a></li><li><a href="#s99">附录：方法与文件</a></li>
  </ol>''', 's00')

# ---------------- 01 physics ----------------
s01 = sec('01', '原理', '螺线管为什么会按质量“分拣”', inner(M, 94, 107), 's01')

# ---------------- 02 evolution ----------------
m03 = inner(M, 109, 117)
m03 = sub(m03, '图 3 是全场追踪的结果', '[[p-cone]] 是全场追踪的结果')
m04 = inner(M, 119, 153)
m04 = sub(m04, '磁屏蔽 B ≲ 10 mT', '第一轮假设：源区屏蔽')
m04 = sub(m04, '<td>位于 z = −0.6 m，正则角动量为零（源区屏蔽）</td>', '<td>位于 z = −0.6 m，正则角动量为零（第一轮假设源区屏蔽，后来证明不需要，见第 06 节）</td>')
m04 = sub(m04, '<td class="num">σ′ = 4 mrad, σ<sub>E</sub> = 1 eV</td><td>径向、方位角散各 4 mrad（rms）</td>', '<td class="num">σ′ = 4 mrad, σ<sub>E</sub> = 1 eV</td><td>径向、方位角散各 4 mrad（rms，假设值；实际设计的源见第 04 节）</td>')
m05 = inner(M, 155, 164)
m06 = inner(M, 166, 171)
m06 = sub(m06, '如果只想做小流强的原理验证，这个构型最简单。', '如果只想做小流强的原理验证，这个构型最简单（第 13 节）。')
s02 = sec('02', '方案概念与第一轮结果', '环形平行束：质量相关透镜', '''  <div class="note"><b>本节是历史背景。</b>这里记录第一轮计算：环形源方案的概念和初步结果。第一轮用的是一个<b>假设的</b>源（角散 4 mrad，源区屏蔽到零场），工作点为 2.74 T、焦面 z = 1.52 m、效率 81%。设计出真实离子源后，工作点变为 2.785 T、焦面 1.365 m、效率约 100%，以第 03 节起的结果为准。</div>
  <h3>方案：环形平行束 → 质量相关透镜 → 中心孔收集（第一轮参数）</h3>
''' + m04 + '''
  <h3>第一轮模拟结果：焦面上的分离</h3>
''' + m05 + '''
  <h3>单束情形：你最初设想的“出射角度差”</h3>
''' + m06, 's02')

# ---------------- 03 layout ----------------
b01 = inner(B, 99, 119)
b01 = sub(b01, '<td>环形离子源（第 13–14 节）</td>', '<td>环形离子源（第 04–05 节）</td>')
b01 = sub(b01, '<td>收集器（见第 06 节）</td>', '<td>收集器（见第 07 节）</td>')
b02 = inner(B, 121, 129)
s03 = sec('03', '整机布局与基准参数', '从离子源到收集器', '''  <p>整机沿 z 轴排列：环形离子源出口在 z = −0.6 m，超导螺线管线圈在 0–1 m，收集器在 z = 1.365 m。下表之后的数值均来自全三维 Boris 追踪，源相空间取自 0.3 mA 自洽引出计算（第 04 节）。</p>
  <div class="kpis">
    <div class="kpi"><b>2.785 T</b><span>工作点（收集器固定时窗口中心），窗口 ±0.18%</span></div>
    <div class="kpi"><b>5.5 → 7.4 mm</b><span>焦面上 Rb-87 最大半径与 Rb-85 最小半径，孔半径取 5.9 mm</span></div>
    <div class="kpi"><b>0.9 m/mA</b><span>未中和流强引起的焦点后移，0.3 mA 为上限</span></div>
    <div class="kpi"><b>≤ 4×10⁻⁷ Pa</b><span>收集器附近铷蒸气分压上限（冷壁 ≤ 250 K）</span></div>
  </div>
  <h3>束流线布局</h3>
''' + b01 + '''
  <h3>束流包络：从源出口到焦面</h3>
''' + b02, 's03')

# ---------------- 04 source optics ----------------
m13 = inner(M, 248, 281)
m13 = sub(m13, '这一节用轴对称静电场求解加离子追踪，设计出能直接产生平行环形束的引出系统，再把源出口的真实相空间分布接入螺线管分离模拟。结论是：源能达到的束流品质优于第一轮的假设（角散约 1.3–2.5 mrad，假设值为 4 mrad），整机在 2.78 T 下 Rb-87 的收集效率接近 100%，丰度 ≥ 99.97%。',
           '这一节用轴对称静电场求解加离子追踪，设计出能直接产生平行环形束的引出系统，再把源出口的真实相空间分布接入螺线管分离模拟。结论是：源能达到的束流品质优于第一轮的假设（角散约 1.3–2.5 mrad，假设值为 4 mrad），整机在 2.78 T 下 Rb-87 的收集效率接近 100%，丰度 ≥ 99.97%。本节的电极是直角的初版，第 05 节给出倒圆角的修订版，出口束流品质不变；束流线计算（第 06 节起）在 2.785 T、孔半径 5.9 mm 下重新确认了这一结果。')
m13 = sub(m13, '<tr><td>磁屏蔽</td><td class="num">源区残余场 ≲ 10 mT</td><td>满足 Busch 条件（见第 07 节）</td></tr>',
           '<tr><td>源区磁场</td><td class="num">发射面处边缘场 16 mT，不屏蔽</td><td>发射面残余场 −10 到 +25 mT 都不影响分离（第 06 节）；可选补偿线圈微调</td></tr>')
m13 = sub(m13, '（已在第 14 节解决：修订后为 3.1 kV/mm）', '（已在第 05 节解决：修订后为 3.1 kV/mm）')
s04 = sec('04', '离子源：引出光学', '环形平行离子源：电极设计与仿真', m13, 's04')

# ---------------- 05 source feasibility ----------------
m14 = inner(M, 283, 357)
m14 = sub(m14, '<p>第 13 节只算了离子光学。', '<p>第 04 节只算了离子光学。')
m14 = sub(m14, '在第 13 节已验证影响可忽略。', '在第 04 节已验证影响可忽略。')
m14 = sub(m14, '<p>第 13 节的模型中电极是直角', '<p>第 04 节的模型中电极是直角')
m14 = sub(m14, '有利于对分离器中的束流做空间电荷中和（第 08 节）。', '有利于对分离器中的束流做空间电荷中和（第 09 节）。')
m14 = sub(m14, '<tr><td>磁屏蔽</td><td class="num">软铁 / 坡莫合金筒</td><td>源区残余场 ≲ 10 mT</td></tr>',
           '<tr><td>补偿线圈（可选）</td><td class="num">半径 80 mm，约 2 kA·匝</td><td>把发射面残余场调零；不加也能工作（第 06 节）</td></tr>')
s05 = sec('05', '离子源：工程可行性', '环形离子源能不能做出来：逐项核算', m14, 's05')

# ---------------- 06 beamline ----------------
b03 = inner(B, 131, 144)
b03 = sub(b03, '第一轮据此要求残余场 ≤ 10 mT。', '第一轮（第 02 节，假设 4 mrad 的源）据此要求残余场 ≤ 10 mT。')
b04 = inner(B, 146, 160)
b04 = sub(b04, '<td class="num">±1 kV（第 13 节）</td>', '<td class="num">±1 kV（第 04 节）</td>')
b05 = inner(B, 162, 184)
b05 = sub(b05, '按图 3 的质谱确认', '按 [[pB]] 的质谱确认')
s06 = sec('06', '束流线', '边缘场、调谐窗口与对准', '''  <h3>源区边缘场和漂移长度</h3>
''' + b03 + '''
  <h3>磁场与高压窗口</h3>
''' + b04 + '''
  <h3>对准误差和静电导向器</h3>
''' + b05, 's06')

# ---------------- 07 collector ----------------
b06 = inner(B, 186, 209)
s07 = sec('07', '收集器与分离性能', '收集器与废束几何', '''  <p>整机分离性能（真实源分布、2.78 T）见第 04 节 [[p-src-sep]]：Rb-87 斑 rms 半径 0.7 mm，Rb-85 环中值半径 9.1 mm。本节按束流线工作点 2.785 T 给出收集器的具体几何。</p>
''' + b06, 's07')

# ---------------- 08 tolerances ----------------
b10 = inner(B, 285, 307)
m07 = inner(M, 173, 184)
m07 = sub(m07, '<div class="note"><b>源区残余磁场（Busch 条件）。</b>', '<div class="note"><b>源区残余磁场（Busch 条件）。</b><span class="upd">已被第 06 节取代</span> ')
m07 = sub(m07, '因此需要一个简单的软铁或坡莫合金屏蔽筒，或者把源后移。</div>', '第一轮据此建议加软铁或坡莫合金屏蔽筒。用真实源重新扫描后（第 06 节），−10 到 +25 mT 都没有影响，所以不需要屏蔽。</div>')
s08 = sec('08', '容差汇总', '误差容差汇总', b10 + '''
  <h3>束流品质要求（第一轮扫描）</h3>
  <p>下面四张图是第一轮对源品质的要求，用的是第 02 节的假设源（4 mrad，2.74 T），用来回答“源至少要多好”。实际设计的源（角散 1.3 / 2.5 mrad）在这些要求之内。</p>
''' + m07, 's08')

# ---------------- 09 space charge + neutralization ----------------
m08 = inner(M, 186, 192)
m15 = M  # alias
m15_diag = R(M, 372, 372)
m15_neut = R(M, 391, 396)
b07 = inner(B, 211, 237)
b07 = sub(b07, '和第 08 节结论一致。', '和第一轮及第 04 节的结论一致。')
s09 = sec('09', '空间电荷与中和', '流强上限、焦点漂移与中和', '''  <p>空间电荷是这台机器最紧的约束。本节先给出第一轮的估计和定性诊断，再给出用真实源分布、从源出口到焦面整段自洽计算的结果，最后讨论中和方案。</p>
  <h3>第一轮估计</h3>
''' + m08 + '''
  <h3>为什么只有 0.3 mA：诊断</h3>
''' + m15_diag + '''
  <h3>全束流线自洽计算</h3>
''' + b07 + '''
  <h3>强磁场中的电子中和</h3>
''' + m15_neut + '''
  <div class="caution"><b>中和度的稳定性要求。</b>上面的焦点灵敏度（约 0.9 m/mA）是对基准设计算的。放大设计（第 11 节）中这一灵敏度应随定标降低，但还没有计算。在得到数值之前，应把“中和度能否稳定在 0.1% 量级”和“中和度能否做到 90% 以上”一起作为中和实验的测量目标。</div>''', 's09')

# ---------------- 10 vacuum ----------------
b08 = inner(B, 239, 268)
m15_cx = R(M, 399, 408)
s10 = sec('10', '真空、电荷交换与铷蒸气', '电荷交换、散射与铷蒸气', b08 + '''
  <h3>放大设计在中和气压下的电荷交换</h3>
''' + m15_cx, 's10')

# ---------------- 11 scaling ----------------
m15_intro = R(M, 362, 369)
m15_A = R(M, 377, 388)
m15_A = sub(m15_A, '<td>基准（第 13–14 节）</td>', '<td>基准（第 03–05 节）</td>')
m15_src = R(M, 411, 416)
m15_src = sub(m15_src, '（图 19 的虚线）', '（[[p-f-cov]] 的虚线）')
m15_col = R(M, 419, 419)
m15_yield = R(M, 422, 432)
m15_tail = re.sub(r'；<a href="https://www.law.cornell.edu[^<]*</a>', '', R(M, 435, 435))
b09 = inner(B, 270, 283)
b09 = sub(b09, '这部分由离子源讨论串负责。', '见本节“离子源放大”。')
s11 = sec('11', '流强放大与产量', '达到电磁分离器的流强与产量', m15_intro + '''
  <h3>两种定标方式</h3>
  <p>保持色散 Δα ∝ ∫B dz/√E 不变时，空间电荷的散焦力与螺线管聚焦力之比约为 K<sub>perv</sub>/(K²R²)，K<sub>perv</sub> ∝ I/E<sup>1.5</sup>，所以可承受的流强满足 I<sub>max</sub> ∝ E<sup>1.5</sup>R²。有两种放大方式，代价不同：</p>
  <div class="tbl"><table>
    <thead><tr><th>方式</th><th>做法</th><th>例子</th><th>无中和上限</th><th>代价</th></tr></thead>
    <tbody>
      <tr><td>A：放大几何</td><td>能量和环半径一起放大，场强保持约 2.8 T，线圈加长、孔径加大</td><td class="num">30 keV, R = 50 mm, 1.7 m × 2.8 T, φ300</td><td class="num">3–4 mA（已模拟）</td><td>磁体更长更粗，场强不变</td></tr>
      <tr><td>C：几何不变</td><td>只提高能量，B ∝ √E，轨迹完全相似</td><td class="num">40 keV, R = 30 mm, 1 m × 5.57 T</td><td class="num">2.4 mA（按 E<sup>1.5</sup> 定标）</td><td>磁体尺寸不变，场强升高；源需 40 kV 引出</td></tr>
    </tbody></table></div>
  <p>两者都符合同一条定标律。方式 A 已用自洽空间电荷模拟验证，方式 C 是严格的相似律推论。两者也可以组合。</p>
  <h3>方式 A：放大几何</h3>
''' + m15_A + '''
  <h3>方式 C：几何不变、提高场强</h3>
''' + b09 + '''
  <h3>离子源放大</h3>
''' + m15_src + '''
  <h3>收集器：一个容易被忽视的问题</h3>
''' + m15_col + '''
  <h3>产量</h3>
''' + m15_yield + '\n' + m15_tail, 's11')

# ---------------- 12 PM ----------------
s12 = sec('12', '永磁版本', '常温永磁螺线管能不能做', inner(M, 194, 206), 's12')

# ---------------- 13 compare + cost ----------------
m10 = inner(M, 208, 224)
m10 = sub(m10, '<td>Larmor 相位 ∝ m<sup>−1/2</sup>，焦面分离约 12 mm（环半径中值）</td>', '<td>Larmor 相位 ∝ m<sup>−1/2</sup>，焦面上 Rb-87 ≤ 5.5 mm、Rb-85 ≥ 7.4 mm</td>')
m10 = sub(m10, '<td>2.74 T、1 m、内径 200 mm 的无液氦超导线圈，冷质量约数百 kg</td>', '<td>2.785 T、1 m、内径 200 mm 的无液氦超导线圈，冷质量约数百 kg</td>')
m10 = sub(m10, '<td>模拟值 99.9%（未计散射）</td>', '<td>模拟值 ≥ 99.97%（真实源；散射和热铷污染预算见第 10 节）</td>')
m10 = sub(m10, '<td>≤ 0.1 mA 已验证；mA 级需要中和，尚未验证</td>', '<td>无中和 0.3 mA（基准）/ 3–4 mA（放大设计）已模拟；更高需要中和，尚未验证</td>')
m11 = inner(M, 226, 235)
m11 = re.sub(r'\n\s*<li><b>永磁替代：</b>.*?</li>', '', m11)
m11 = sub(m11, '<li><b>产量：</b>1 mA Rb⁺ 总流（含两种同位素）× 27.8% × 81% ≈ 1.4×10¹⁵ 个 Rb-87/s ≈ <b>17 mg/天</b>，按全年连续运行约 6 g。一个铷原子钟气室只需毫克量级的 Rb-87，因此即使只有 0.1 mA，每天 1.7 mg 也足以满足原子钟的用量。</li>',
           '<li><b>产量：</b>基准设计在无中和上限 0.3 mA 下约 <b>6 mg/天</b> Rb-87（0.3 mA × 27.8% × 收集效率约 100%），全年连续运行约 2 g。一个铷原子钟气室只需毫克量级的 Rb-87，因此基准设计已足以满足原子钟的用量；更高产量见第 11 节。</li>')
m11 = sub(m11, '<li><b>磁体：</b>', '<li><b>离子源：</b>源本体加电源粗估 30–60 万元人民币（第 05 节）。</li>\n    <li><b>磁体：</b>')
s13 = sec('12', '对比与成本', '与传统二极磁铁电磁分离器比较', m10 + '''
  <h3>成本与产量估算</h3>
''' + m11, 's13')

# ---------------- 14 next steps ----------------
s14 = sec('13', '风险与下一步', '建议的验证路线', '''  <ul class="tight">
    <li><b>中和实验（最关键）。</b>在现有 1–3 T 螺线管上，用 0.1–10 mA 的 Rb⁺ 束，在 10⁻⁶–10⁻⁴ torr 范围内扫描氩气压，并试验电子注入（第 09 节），测量焦点位置、焦斑和丰度随气压的变化，以及中和度的稳定性。</li>
    <li><b>单束原理验证。</b>用第 02 节的单束构型，以 0.1–1 µA 的 Rb⁺ 测量焦面分布和丰度，验证 Larmor 色散的绝对值。</li>
    <li><b>离子源样件。</b>做一个小样件（直条形或小半径环），实测 Re/Ir 表面的电离效率、发射度和多孔钨透过率；同时用 IBSimu 或 COMSOL 做三维校核，评估加工误差造成的方位不均匀（第 05 节）。</li>
    <li><b>收集器三维校核。</b>用静电求解器对带孔的减速口袋和刀口孔嘴做一次三维计算，确认口袋的场不会穿过孔口偏转 Rb-85（第 07 节）。</li>
    <li><b>线圈误差场。</b>采购时测量线圈的实际误差场（制造公差带来的非轴对称谐波），用本项目代码算一次。</li>
    <li><b>多误差蒙特卡洛。</b>目前容差都是逐项单独加入，还没有做多项误差同时加入的计算（第 08 节）。</li>
    <li><b>放大设计的焦点灵敏度。</b>计算放大设计 A 中焦点随未中和流强的漂移，确定其中和度稳定性要求（第 09 节）。</li>
  </ul>
  <div class="caution"><b>仍未覆盖的风险。</b>（1）材料参数（功函数、Rb⁺ 脱附能、多孔钨透过率、电荷交换与散射截面）都取自通用文献量级，相关结果只有数量级意义。（2）所有静电计算都是二维轴对称的。（3）强磁场中离子–电子等离子体可能激发二极管或离子回旋类不稳定性，目前无法可靠计算。</div>''', 's14')

# ---------------- appendix ----------------
m_app = inner(M, 438, 448)
b_app = inner(B, 319, 324)
b_app = sub(b_app, '磁场、追踪方法同主报告（250×4 电流环，椭圆积分，Boris 推进器，ω<sub>c</sub>Δt = 0.02）。', '')
b_app = sub(b_app, '源相空间重新运行了第 13 节的引出计算（0.3 mA 自洽，6000 个粒子；第 14 节的圆角电极修订不改变出口束流品质）。', '<b>束流线计算（第 03、06–10 节）：</b>源相空间重新运行了第 04 节的引出计算（0.3 mA 自洽，6000 个粒子；第 05 节的圆角电极修订不改变出口束流品质）。')
app = '''<section id="s99">
  <div class="secno">附录 · 方法与文件</div>
  <h2>模拟方法</h2>
''' + m_app + '\n' + b_app + '''
  <p>其余代码和数据：<span class="k">solenoid-separator/sim/</span>（第一轮：分离、容差、空间电荷）、<span class="k">solenoid-separator/source/</span>（离子源光学与可行性，第 04–05 节）、<span class="k">solenoid-separator/highcurrent/</span>（流强放大与电荷交换蒙特卡罗，第 09–11 节）。</p>
</section>
'''

body = header + '\n' + '\n'.join([summary, s01, s02, s03, s04, s05, s06, s07, s08, s09, s10, s11, s13, s14, app]) + '</div>\n'

# ---------------- figure renumbering ----------------
fig_re = re.compile(r'(<figure><div class="plot[^"]*" id="([^"]+)"></div><figcaption><span class="tag">)图 \d+(</span>)')
num = {}
def repl(m):
    num[m.group(2)] = len(num) + 1
    return f'{m.group(1)}图 {num[m.group(2)]}{m.group(3)}'
body = fig_re.sub(repl, body)
for k, v in num.items():
    body = body.replace(f'[[{k}]]', f'图 {v}')
left = re.findall(r'\[\[[^\]]+\]\]', body)
if left: sys.exit('unresolved refs ' + str(left))
if re.search(r'第 1[3-5] 节|第 13–14 节', body): print('WARN stale ref', re.findall(r'.{20}第 1[3-5] 节', body))

# ---------------- scripts ----------------
ms = next(i for i, l in enumerate(M) if l.startswith('<script src="https://cdnjs'))
me = max(i for i, l in enumerate(M) if l.startswith('</script>'))
mscript = '\n'.join(M[ms:me+1])
mscript = mscript.replace('第 13 节', '第 04 节')
mscript = re.sub(r'  // cone\n.*?(?=  // envelope)', '', mscript, flags=re.S)

FIX = """<script>
/* Plotly 的 clip-path 用 url(#id) 引用；页面带 <base> 时部分浏览器（Safari、Firefox）解析到别的地址，SVG 数据层被整体裁掉。改成带当前页地址的绝对引用。 */
(function(){
  if(!window.Plotly||location.protocol==='about:')return;
  const page=location.href.split('#')[0];
  function fixClip(gd){gd.querySelectorAll('[clip-path]').forEach(e=>{const m=(e.getAttribute('clip-path')||'').match(/^url\\(['"]?#([^'")]+)['"]?\\)$/);if(m)e.setAttribute('clip-path','url("'+page+'#'+m[1]+'")')});}
  const react=Plotly.react;
  Plotly.react=function(){return react.apply(Plotly,arguments).then(gd=>{fixClip(gd);if(!gd.__fixClip){gd.__fixClip=1;gd.on('plotly_afterplot',()=>fixClip(gd));}return gd;});};
})();
</script>"""

if 'k-T' in mscript: print('note: k-T referenced in script')
bs = next(i for i, l in enumerate(B) if l.startswith('<script>'))
be = max(i for i, l in enumerate(B) if l.startswith('</script>'))
bscript = '\n'.join(B[bs+1:be])
bscript = re.sub(r'^const D=', 'const DB=', bscript, flags=re.M)
bscript = re.sub(r'(?<![\w.$])D(?=[.\[])', 'DB', bscript)
bscript = '\n'.join(l for l in bscript.split('\n') if not l.startswith('const css='))
bscript = re.sub(r'\bfunction plot\(', 'function svgPlot(', bscript)
bscript = re.sub(r'(?<![\w.])plot\(\'', "svgPlot('", bscript)
bscript = re.sub(r'\bfunction draw\(\)', 'function drawBL()', bscript)
bscript = re.sub(r'(?<![\w.])draw\b(?!\w)', 'drawBL', bscript).replace('function drawBL()', 'function drawBL()')
bscript = bscript.replace("matchMedia('(prefers-color-scheme: dark)').addEventListener('change',drawBL);","try{matchMedia('(prefers-color-scheme: dark)').addEventListener('change',drawBL)}catch(e){}")
bscript = '<script>\n' + bscript + '\n</script>'

out = head + '\n' + body + '\n' + mscript + '\n' + bscript + '\n'
open('report.html', 'w').write(out)
print('ok', len(out), 'figs', len(num))
