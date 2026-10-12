"""Pictograms for the modes, drawn as SVG in the current theme's colours."""
from PySide6.QtCore import QByteArray, QRectF, Qt
from PySide6.QtGui import QPainter, QPixmap
from PySide6.QtSvg import QSvgRenderer

# Each pictogram uses {a} (accent), {b} (secondary lines) and {f} (soft accent fill) on a 48 x 48 grid.
_HEADS = ''.join(f'<circle cx="{x}" cy="{{y}}" r="2.6" fill="{{a}}"/>' for x in range(6, 46, 6))

PICTOGRAMS = {
    # Transplant: a protein lifted from one system into another.
    'charmprot': '''
<rect x="3" y="15" width="17" height="24" rx="3" fill="none" stroke="{b}" stroke-width="1.8"/>
<rect x="28" y="15" width="17" height="24" rx="3" fill="none" stroke="{b}" stroke-width="1.8"/>
<rect x="8" y="21" width="7" height="12" rx="3.5" fill="{f}" stroke="{a}" stroke-width="1.8"/>
<rect x="33" y="21" width="7" height="12" rx="3.5" fill="{a}"/>
<path d="M12 12 C 16 3, 32 3, 36 12" fill="none" stroke="{a}" stroke-width="2" stroke-linecap="round"/>
<path d="M32.5 10.5 L36 12.5 L37 8.6" fill="none" stroke="{a}" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"/>''',
    # Protein insertion: a helix bundle entering a bilayer.
    'protein': '''
<rect x="2" y="22" width="44" height="14" rx="2" fill="{f}"/>
<line x1="2" y1="22" x2="46" y2="22" stroke="{b}" stroke-width="1.6" stroke-dasharray="2 3"/>
<line x1="2" y1="36" x2="46" y2="36" stroke="{b}" stroke-width="1.6" stroke-dasharray="2 3"/>
<rect x="15" y="8" width="5" height="34" rx="2.5" fill="{a}"/>
<rect x="21.5" y="5" width="5" height="34" rx="2.5" fill="{a}" opacity="0.8"/>
<rect x="28" y="9" width="5" height="34" rx="2.5" fill="{a}" opacity="0.6"/>''',
    # Ligand swap: one ring out, another in.
    'lig': '''
<path d="M12 9 L19 13 L19 21 L12 25 L5 21 L5 13 Z" fill="none" stroke="{b}" stroke-width="1.8" stroke-linejoin="round"/>
<path d="M36 23 L43 27 L43 35 L36 39 L29 35 L29 27 Z" fill="{f}" stroke="{a}" stroke-width="2" stroke-linejoin="round"/>
<circle cx="36" cy="31" r="3.2" fill="none" stroke="{a}" stroke-width="1.6"/>
<line x1="43" y1="27" x2="46" y2="23" stroke="{a}" stroke-width="2" stroke-linecap="round"/>
<path d="M22 12 C 31 12, 34 15, 34 20" fill="none" stroke="{a}" stroke-width="1.8" stroke-linecap="round"/>
<path d="M31.5 18 L34 21 L36.5 18" fill="none" stroke="{a}" stroke-width="1.8" stroke-linecap="round" stroke-linejoin="round"/>
<path d="M26 36 C 17 36, 14 33, 14 28" fill="none" stroke="{b}" stroke-width="1.8" stroke-linecap="round"/>
<path d="M11.5 30 L14 27 L16.5 30" fill="none" stroke="{b}" stroke-width="1.8" stroke-linecap="round" stroke-linejoin="round"/>''',
    # Cholesterol: the fused steroid rings with the hydroxyl and the tail.
    'chl': '''
<path d="M5 30 L10 27 L15 30 L15 36 L10 39 L5 36 Z" fill="{f}" stroke="{a}" stroke-width="1.8" stroke-linejoin="round"/>
<path d="M15 30 L20 27 L25 30 L25 36 L20 39 L15 36" fill="none" stroke="{a}" stroke-width="1.8" stroke-linejoin="round"/>
<path d="M20 27 L20 21 L25 18 L30 21 L30 27 L25 30" fill="none" stroke="{a}" stroke-width="1.8" stroke-linejoin="round"/>
<path d="M30 21 L35 18 L38.5 22 L35.5 26.5 L30 27" fill="none" stroke="{a}" stroke-width="1.8" stroke-linejoin="round"/>
<path d="M35 18 L37 12 L42 10 L44 5" fill="none" stroke="{b}" stroke-width="1.8" stroke-linecap="round" stroke-linejoin="round"/>
<line x1="5" y1="36" x2="2.5" y2="41" stroke="{b}" stroke-width="1.8" stroke-linecap="round"/>
<circle cx="2.8" cy="42.5" r="2" fill="{b}"/>''',
    # Add binder: a receptor in the membrane, a ligand above it.
    'addbinder': '''
<rect x="2" y="30" width="44" height="12" rx="2" fill="{f}"/>
<line x1="2" y1="30" x2="46" y2="30" stroke="{b}" stroke-width="1.6" stroke-dasharray="2 3"/>
<line x1="2" y1="42" x2="46" y2="42" stroke="{b}" stroke-width="1.6" stroke-dasharray="2 3"/>
<rect x="17" y="18" width="14" height="28" rx="6" fill="{a}"/>
<line x1="24" y1="7.5" x2="24" y2="16" stroke="{b}" stroke-width="1.6" stroke-dasharray="1.5 2.5" stroke-linecap="round"/>
<path d="M24 1 L29 4 L29 9.5 L24 12.5 L19 9.5 L19 4 Z" fill="none" stroke="{a}" stroke-width="2" stroke-linejoin="round"/>''',
    # POPC: a bilayer, head groups out and tails in.
    'popc': '''
<g stroke="{b}" stroke-width="1.5" stroke-linecap="round" fill="none">
''' + ''.join(f'<path d="M{x - 0.8} 11 L{x - 0.8} 22 M{x + 0.8} 11 L{x + 0.8} 16 L{x + 2} 18 L{x + 0.8} 22"/>'
              f'<path d="M{x - 0.8} 37 L{x - 0.8} 26 M{x + 0.8} 37 L{x + 0.8} 32 L{x + 2} 30 L{x + 0.8} 26"/>'
              for x in range(6, 46, 6)) + '''
</g>
''' + _HEADS.replace('{y}', '8') + _HEADS.replace('{y}', '40'),
    # Membrane: a bilayer of mixed lipids with cholesterol rings among them.
    'membrane': '''
<rect x="2" y="12" width="44" height="24" rx="2" fill="{f}"/>
''' + ''.join(f'<circle cx="{x}" cy="9" r="2.6" fill="{"{a}" if n % 2 else "{b}"}"/>'
              f'<circle cx="{x}" cy="39" r="2.6" fill="{"{b}" if n % 2 else "{a}"}"/>'
              for n, x in enumerate(range(6, 46, 6))) + '''
<path d="M17 15 L21 17 L21 21 L17 23 L13 21 L13 17 Z" fill="none" stroke="{a}" stroke-width="1.6"/>
<path d="M33 25 L37 27 L37 31 L33 33 L29 31 L29 27 Z" fill="none" stroke="{a}" stroke-width="1.6"/>''',
    # Salt: cations and anions in water.
    'salt': '''
<circle cx="13" cy="14" r="7" fill="{a}"/>
<path d="M9.5 14 H16.5 M13 10.5 V17.5" stroke="{f}" stroke-width="2" stroke-linecap="round"/>
<circle cx="34" cy="20" r="9" fill="none" stroke="{b}" stroke-width="2"/>
<path d="M29.5 20 H38.5" stroke="{b}" stroke-width="2" stroke-linecap="round"/>
<circle cx="19" cy="35" r="9" fill="none" stroke="{b}" stroke-width="2"/>
<path d="M14.5 35 H23.5" stroke="{b}" stroke-width="2" stroke-linecap="round"/>
<circle cx="38" cy="38" r="6" fill="{a}"/>
<path d="M35 38 H41 M38 35 V41" stroke="{f}" stroke-width="2" stroke-linecap="round"/>''',
    # Protonation: an imidazole ring picking up a proton.
    'protonate': '''
<path d="M18 16 L28 20 L27 31 L16 34 L10 25 Z" fill="{f}" stroke="{a}" stroke-width="2" stroke-linejoin="round"/>
<line x1="10" y1="25" x2="3" y2="25" stroke="{b}" stroke-width="1.8" stroke-linecap="round"/>
<path d="M32 13 C 36 9, 40 9, 41 12" fill="none" stroke="{b}" stroke-width="1.6" stroke-dasharray="2 2"/>
<line x1="28" y1="20" x2="31" y2="16" stroke="{a}" stroke-width="2" stroke-linecap="round"/>
<circle cx="40" cy="12" r="6" fill="{a}"/>
<path d="M37 12 H43 M40 9 V15" stroke="{f}" stroke-width="2" stroke-linecap="round"/>''',
    # Mutation: one residue of a chain swapped for another.
    'mutate': '''
<path d="M3 30 C 9 22, 13 38, 19 30 S 29 22, 35 30 S 43 38, 46 30" fill="none" stroke="{b}" stroke-width="2"
      stroke-linecap="round"/>
<line x1="19" y1="30" x2="19" y2="18" stroke="{b}" stroke-width="1.8"/>
<circle cx="19" cy="15" r="4" fill="none" stroke="{b}" stroke-width="1.8"/>
<line x1="35" y1="30" x2="35" y2="17" stroke="{a}" stroke-width="2"/>
<path d="M35 6 L40 9 L40 15 L35 18 L30 15 L30 9 Z" fill="{f}" stroke="{a}" stroke-width="2" stroke-linejoin="round"/>
<path d="M23 12 C 26 9, 27 9, 29 10" fill="none" stroke="{a}" stroke-width="1.6" stroke-linecap="round"/>''',
    # HMR: mass moves from a heavy atom to its hydrogens.
    'hmr': '''
<circle cx="24" cy="24" r="9" fill="{f}" stroke="{a}" stroke-width="2"/>
<line x1="24" y1="15" x2="24" y2="9" stroke="{b}" stroke-width="1.8"/>
<line x1="16" y1="29" x2="10" y2="33" stroke="{b}" stroke-width="1.8"/>
<line x1="32" y1="29" x2="38" y2="33" stroke="{b}" stroke-width="1.8"/>
<circle cx="24" cy="6" r="4.2" fill="{a}"/>
<circle cx="8" cy="35" r="4.2" fill="{a}"/>
<circle cx="40" cy="35" r="4.2" fill="{a}"/>''',
    # Resize: a box made taller, the new water above and below.
    'resize': '''
<rect x="10" y="10" width="28" height="28" rx="2" fill="none" stroke="{b}" stroke-width="1.8"/>
<rect x="10" y="20" width="28" height="8" fill="{a}" opacity="0.85"/>
<rect x="10" y="3" width="28" height="5" rx="1" fill="{f}" stroke="{a}" stroke-width="1.4" stroke-dasharray="2 2"/>
<rect x="10" y="40" width="28" height="5" rx="1" fill="{f}" stroke="{a}" stroke-width="1.4" stroke-dasharray="2 2"/>
<path d="M43 12 V36 M40.5 15 L43 12 L45.5 15 M40.5 33 L43 36 L45.5 33" fill="none" stroke="{a}" stroke-width="1.8"
      stroke-linecap="round" stroke-linejoin="round"/>''',
    # Flood: many ligand copies scattered in the water above a membrane.
    'flood': '''
<rect x="2" y="34" width="44" height="10" rx="2" fill="{f}"/>
<line x1="2" y1="34" x2="46" y2="34" stroke="{b}" stroke-width="1.6" stroke-dasharray="2 3"/>
''' + ''.join(f'<path d="M{x} {y - 4} L{x + 3.5} {y - 2} L{x + 3.5} {y + 2} L{x} {y + 4} L{x - 3.5} {y + 2} '
              f'L{x - 3.5} {y - 2} Z" fill="none" stroke="{{a}}" stroke-width="1.7" stroke-linejoin="round"/>'
              for x, y in ((8, 8), (24, 6), (39, 10), (15, 21), (32, 23), (42, 26))),
    # Oligomer: two receptors side by side in one membrane.
    'oligomer': '''
<rect x="2" y="18" width="44" height="14" rx="2" fill="{f}"/>
<line x1="2" y1="18" x2="46" y2="18" stroke="{b}" stroke-width="1.6" stroke-dasharray="2 3"/>
<line x1="2" y1="32" x2="46" y2="32" stroke="{b}" stroke-width="1.6" stroke-dasharray="2 3"/>
<rect x="7" y="9" width="13" height="32" rx="6" fill="{a}"/>
<rect x="28" y="9" width="13" height="32" rx="6" fill="none" stroke="{a}" stroke-width="2"/>
<path d="M21.5 25 H26.5" stroke="{b}" stroke-width="1.6" stroke-dasharray="1.5 2"/>''',
    # Minimization: down an energy landscape into its minimum.
    'minimize': '''
<line x1="5" y1="4" x2="5" y2="43" stroke="{b}" stroke-width="1.6" stroke-linecap="round"/>
<line x1="5" y1="43" x2="45" y2="43" stroke="{b}" stroke-width="1.6" stroke-linecap="round"/>
<path d="M8 8 C 14 8, 15 22, 20 22 C 24 22, 24 15, 28 15 C 33 15, 32 37, 37 37 C 41 37, 42 30, 45 26"
      fill="none" stroke="{a}" stroke-width="2.2" stroke-linecap="round"/>
<circle cx="37" cy="33" r="3.4" fill="{a}"/>
<path d="M12 6 L16 12" stroke="{b}" stroke-width="1.4" stroke-dasharray="1.5 2" stroke-linecap="round"/>''',
}


def pictogram(mode, accent, secondary, fill, size=40, ratio=2.0):
    """A crisp QPixmap of the mode's pictogram in the given colours."""
    body = PICTOGRAMS[mode].replace('{a}', accent).replace('{b}', secondary).replace('{f}', fill)
    svg = f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 48 48">{body}</svg>'
    renderer = QSvgRenderer(QByteArray(svg.encode()))
    pixmap = QPixmap(int(size * ratio), int(size * ratio))
    pixmap.fill(Qt.GlobalColor.transparent)
    painter = QPainter(pixmap)
    painter.setRenderHint(QPainter.RenderHint.Antialiasing)
    renderer.render(painter, QRectF(0, 0, size * ratio, size * ratio))
    painter.end()
    pixmap.setDevicePixelRatio(ratio)
    return pixmap
