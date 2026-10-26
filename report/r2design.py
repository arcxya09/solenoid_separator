# 第二轮评审核验结果写入设计报告：design_body_B.html -> design_body_B2.html
# 用法：python3 r2design.py rated   （额定 5 mA / 极限 6 mA）
#       python3 r2design.py six     （保持 6 mA 设计点）
import sys
mode = sys.argv[1]
s = open('design_body_B.html').read()
def R(old, new, n=1):
    global s
    if s.count(old) != n: sys.exit(f'count {s.count(old)}: ' + old[:80])
    s = s.replace(old, new)

if mode == 'rated':
    R('设计流强 6 mA，不需要空间电荷中和。', '额定流强 5 mA、极限 6 mA，不需要空间电荷中和。')
    R('<div class="kpi"><b>6 mA</b><span>设计束流（不中和），Rb⁺ 30 keV</span></div>', '<div class="kpi"><b>5 / 6 mA</b><span>额定 / 极限束流（不中和），Rb⁺ 30 keV</span></div>')
    R('<div class="kpi"><b>0.13 g/天</b><span>Rb-87 产量（6 mA 连续运行，约 47 g/年）</span></div>', '<div class="kpi"><b>0.11–0.13 g/天</b><span>Rb-87 产量（5–6 mA 连续运行，约 40–47 g/年）</span></div>')
    R('<td class="num">6 mA，不中和（源流强稳定度 ±5%）</td>', '<td class="num">额定 5 mA，极限 6 mA，不中和（源流强稳定度 ±5%）</td>')
    R('<td>设计流强</td>', '<td>流强</td>')
    R('<tr><td class="num">5 mA</td><td class="num">1.912 T</td><td class="num">2.415 m</td><td class="num">21 mm</td><td class="num">100%，z 窗口 2.36–2.47 m</td></tr>',
      '<tr><td class="num"><b>5 mA（额定）</b></td><td class="num"><b>1.912 T</b></td><td class="num"><b>2.415 m</b></td><td class="num"><b>21 mm</b></td><td class="num"><b>100%，z 窗口 2.36–2.47 m</b></td></tr>')
    R('<tr><td class="num"><b>6 mA（设计点）</b></td><td class="num"><b>1.912 T</b></td><td class="num"><b>2.415 m</b></td><td class="num"><b>21 mm</b></td><td class="num"><b>100%，z 窗口约 40 mm</b></td></tr>',
      '<tr><td class="num">6 mA（极限）</td><td class="num">1.912 T</td><td class="num">2.415 m</td><td class="num">21 mm</td><td class="num">100%，z 窗口约 40 mm</td></tr>')
    concl = '<p><b>运行流强：</b>在设计公差下，额定 5 mA 运行的效率 ≥ 99%；6 mA 运行时 94% 的样本效率 ≥ 95%（中位 98.6%），丰度判据全部样本可满足。收集器 x/y/z 平移台为必配（固定收集器合格率仅 6–25%）。</p>'
else:
    concl = '<p><b>运行流强：</b>在设计公差下，6 mA 运行时 94% 的样本效率 ≥ 95%（中位 98.6%），丰度判据全部样本可满足；5 mA 运行时效率 ≥ 99%。收集器 x/y/z 平移台为必配（固定收集器合格率仅 6–25%）。</p>'

# 源区磁场补偿
R('<tr><td>源区磁场</td><td class="num">发射面 B<sub>z</sub> = 14.5 mT，源出口 23 mT（已计入追踪）</td></tr>',
  '<tr><td>源区磁场</td><td class="num">未补偿时发射面 14.5 mT、源出口 23 mT；三环补偿后源区 |B<sub>z</sub>|、|B<sub>r</sub>| ≤ 0.6 mT</td></tr>')
R('<li><b>源区磁场：</b>正电压单透镜处的残余场保持 ≤ 5 mT（源区磁屏蔽或补偿线圈），防止单透镜在弱场中俘获电子；分离光学对发射面残余场不敏感。</li>',
  '<li><b>源区磁场补偿：</b>正电压单透镜在弱场中会俘获电子，源区残余场须 ≤ 5 mT。采用三组半径 0.25 m 的同轴铜线圈：z = −1.40 / −1.18 / −0.96 m，安匝 −2959 / −2543 / −8142 A·匝。在 r ≤ 0.12 m、z = −1.31 … −1.05 m 区域内 |B<sub>z</sub>| ≤ 0.57 mT、|B<sub>r</sub>| ≤ 0.53 mT；线圈电流稳定度 5×10⁻³ 即可。线圈半径须明显大于束环（≥ 250 mm）并沿轴向分布，单个线圈无法同时压低 B<sub>z</sub> 和 B<sub>r</sub>。替代方案为 3–5 mm 厚软铁屏蔽筒（装配后需实测场图并微调 B₀）。</li>')
# 导向器
R('<td>长 87 mm；8 电极圆柱导向器，或板距 ≥ 280 mm 平板（±10 mrad 需 6.9 kV/m，280 mm 板距为 ±0.97 kV）</td>',
  '<td>长 87 mm，8 电极圆柱导向器（±10 mrad 需 6.9 kV/m）</td>')
R('<td>同 ST1，位于 60–90 mT 边缘场</td><td>与 ST1 一起把环对到磁轴（位置 + 角度）</td>',
  '<td>同 ST1，孔径 140 mm，位于 60–90 mT 边缘场；8 个电极上可叠加 ±100 V 的 cos 2θ / sin 2θ 四极偏压</td><td>与 ST1 一起把环对到磁轴（位置 + 角度）；四极偏压现场微调线圈四极误差</td>')
# 收集器与废束靶
R('<td class="num">石墨，外半径 185 mm，冷媒 ≤ 210 K</td>', '<td class="num">深槽百叶石墨靶（槽深 25 mm、节距 5 mm），外半径 185 mm，冷媒 ≤ 210 K</td>')
R('外锥半角 35°、长 25 mm，与废束靶面相接；内锥斜率 0.30。', '外锥半角 25°，与废束靶面相接；内锥斜率 0.30。外锥半角取 25°（而非 35°）使打在锥外壁上的 Rb-85 从 5.8% 降到 0.23%。')
R('<li><b>功率：</b>Rb-85 束流 4.33 mA × 30 keV = 130 W。94% 落在靶面 r = 40–159 mm（中位 63 mm），共 122 W，峰值 1.4 W/cm²（r = 40–50 mm）；6% 落在锥外壁 r = 28–40 mm，7.5 W，峰值 0.34 W/cm²。</li>\n    <li><b>结构：</b>石墨靶厚 10 mm（温升约 1.4 K），贴在制冷铜背板上，外半径 185 mm；冷媒 ≤ 210 K，靶面 ≤ 230 K。石墨以插层化合物形式固定铷。</li>',
  '<li><b>结构：</b>深槽百叶靶。30 片同心锥面薄片，片顶在同一平面，半径 39–184 mm，节距 5 mm，槽深 25 mm；每片沿当地离子方向外倾（与轴夹 11.4°–21°），刃口 ≤ 0.1 mm，向根部渐厚到约 1 mm。槽底为平板，兼作散热板，贴在制冷铜背板上；外壁 r = 185 mm。Rb-85 打在槽底，溅射原子须经多次冷壁碰撞才能逃出槽口。</li>\n    <li><b>功率：</b>Rb-85 束流 4.33 mA × 30 keV = 130 W（6 mA）。槽底 126 W，入射角中位 14°，峰值 1.03 W/cm²；刃尖约 2 W，片面约 0.5 W，片的温升 &lt; 1 K；锥外壁约 0.3 W。冷媒 ≤ 210 K，所有表面 ≤ 230 K。</li>')
i = s.index('  <h3>冷表面</h3>'); j = s.index('</section>', i)
s = s[:i] + '''  <h3>溅射中性原子</h3>
  <p>溅射中性原子的蒙特卡罗（机型 B 6 mA 自洽束流，每种工况 10⁷ 个原子，斜入射产额按 cos<sup>−1.7</sup> 加权，原子在各表面按粘附系数 s 留下或漫反射）给出每个入射 Rb-85 离子对应的进孔原子数 F：</p>
  <div class="tbl"><table><thead><tr><th>废束靶</th><th>s = 0.85</th><th>s = 0.9</th><th>s = 0.95</th><th>s = 0.99</th></tr></thead><tbody>
    <tr><td>深槽百叶靶 + 外锥 25°（定型）</td><td class="num">3.5×10⁻⁵</td><td class="num">2.2×10⁻⁵</td><td class="num">1.3×10⁻⁵</td><td class="num">1.2×10⁻⁶</td></tr>
    <tr><td>平面石墨靶 + 外锥 35°（对比）</td><td class="num">4.5×10⁻⁴</td><td class="num">3.0×10⁻⁴</td><td class="num">1.4×10⁻⁴</td><td class="num">2.8×10⁻⁵</td></tr>
  </tbody></table></div>
  <p>丰度损失 ≈ 2.6 × Y × F（Y 为 30 keV 自溅射产额）。定型结构下，Y = 2 时为 1.8×10⁻⁴ – 6×10⁻⁶，Y = 10 时为 9×10⁻⁴ – 3×10⁻⁵，在粘附系数低到 0.85（超热溅射原子的动力学反射、表面被金属铷覆盖）时仍满足 99.9%。所有能看到孔口的表面（锥、靶、衬管末端、E2）保持 ≤ 230 K。</p>
''' + s[j:]
R('<td>废束靶溅射</td><td>冷表面 s ≥ 0.99</td><td class="num">~1×10⁻⁴</td>', '<td>废束靶溅射</td><td>深槽百叶靶，s = 0.85–0.99</td><td class="num">≤ 1.8×10⁻⁴（Y = 2）</td>')

# 容差
i = s.index('<section id="d07">'); j = s.index('  <h3>对中步骤</h3>', i)
s = s[:i] + '''<section id="d07">
  <div class="secno">07 · 容差</div>
  <h2>误差容差</h2>
  <p>收集器 x/y/z 可调（z ±50 mm），孔半径 21–22 mm。四极误差 b₂、a₂ 按束流半径 R₀ = 104 mm 处的相对场定义。</p>
  <div class="tbl"><table>
    <thead><tr><th>误差项</th><th>设计公差（1σ）</th><th>宽松公差组（1σ，对比）</th><th>单项影响（6 mA）</th></tr></thead>
    <tbody>
      <tr><td>线圈偏移</td><td class="num">0.5 mm</td><td class="num">0.5 mm</td><td class="num">2σ / 4σ：98.9% / 92.1%</td></tr>
      <tr><td>线圈倾斜</td><td class="num">0.25 mrad</td><td class="num">0.5 mrad</td><td class="num">2σ / 4σ：98.6% / 95.1%（0.5 mrad 时）</td></tr>
      <tr><td>线圈四极场 b₂、a₂（@104 mm）</td><td class="num">2.5×10⁻⁴（线圈指标上限 5×10⁻⁴）</td><td class="num">5×10⁻⁴</td><td class="num">5×10⁻⁴ 时 99.5%，1×10⁻³ 时 97.8%</td></tr>
      <tr><td>源偏移</td><td class="num">0.5 mm</td><td class="num">0.5 mm</td><td class="num">2σ / 4σ：98.0% / 96.1%</td></tr>
      <tr><td>源倾斜</td><td class="num">0.75 mrad</td><td class="num">1.5 mrad</td><td class="num">1 mrad 的 2σ / 4σ：95.5% / 89.7%</td></tr>
      <tr><td>源环半径</td><td class="num">0.25 mm</td><td class="num">0.25 mm</td><td class="num">≥ 99.9%</td></tr>
      <tr><td>ΔB/B、ΔE/E</td><td class="num">3×10⁻⁴</td><td class="num">3×10⁻⁴</td><td class="num">≥ 99.9%</td></tr>
      <tr><td>单透镜电压</td><td class="num">30 V</td><td class="num">100 V</td><td class="num">—</td></tr>
      <tr><td>收集器横向 / 纵向</td><td class="num">0.3 mm / 5 mm</td><td class="num">0.3 mm / 5 mm</td><td class="num">—</td></tr>
    </tbody>
  </table></div>
  <p>单项误差下丰度都 ≥ 99.9%，误差只影响收集效率。线圈四极场的采购指标为上限 ≤ 5×10⁻⁴（104 mm 处），统计上对应 σ = 2.5×10⁻⁴；按束流半径规定（按 30 mm 参考半径规定会放松约 3.5 倍）。ST2 上的 ±100 V 四极偏压可在运行流强下现场微调残余四极误差（收益 0.2–0.7%），但不能代替线圈公差。</p>
  <h3>多误差蒙特卡罗</h3>
  <p>所有误差按上表同时随机抽取（每组 32 个样本，自洽空间电荷，Poisson 迭代 6 次）。另一组独立计算（20 个样本，b₂、a₂ σ = 2.5×10⁻⁴、源倾斜 σ = 1 mrad）给出 5 mA 合格率 100%（平均效率 99.98%），6 mA 平均效率 96.6%。合格判据：Rb-87 效率 ≥ 99% 且丰度 ≥ 99.9%。</p>
  <div class="tbl"><table>
    <thead><tr><th>情形（收集器可调）</th><th>合格率</th><th>效率 ≥ 95% 的比例</th><th>效率中位数</th><th>效率 P10</th></tr></thead>
    <tbody>
      <tr><td>设计公差，6 mA</td><td class="num">50%</td><td class="num">94%</td><td class="num">98.6%</td><td class="num">96.3%</td></tr>
      <tr><td>宽松公差组，5 mA</td><td class="num">69%</td><td class="num">94%</td><td class="num">99.6%</td><td class="num">96.4%</td></tr>
      <tr><td>宽松公差组，6 mA</td><td class="num">28%</td><td class="num">53%</td><td class="num">95.6%</td><td class="num">90%</td></tr>
    </tbody>
  </table></div>
  ''' + concl + '''
  <p>设计公差下影响最大的是源与线圈的横向对中（源偏移、线圈偏移、源倾斜），由导向器 ST1 / ST2 和 B₀ 微调在线校正；上表未计入在线重调，结果偏保守。</p>
''' + s[j:]

# 电源
R('<tr><td>导向器 ST1 / ST2</td><td class="num">±1 kV 量级</td><td>—</td></tr>',
  '<tr><td>导向器 ST1 / ST2</td><td class="num">8 电极，偶极 ±1 kV 量级；ST2 叠加四极 ±100 V</td><td>—</td></tr>\n      <tr><td>源区补偿线圈（3 组）</td><td class="num">3–8 kA·匝，常规铜线圈</td><td>稳定度 5×10⁻³</td></tr>')
open('design_body_B2.html', 'w').write(s)
print('ok', mode)
