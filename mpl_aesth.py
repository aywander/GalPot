import matplotlib.pyplot as plt
from cycler import cycler

global c

def adjust_rcParams(style='seaborn', use_kpfonts=False, dark_mode=False):

    plt.style.use(style)

    if dark_mode:
        fc = 'white'
        fc_i = 'black'
        fc_n = (1, 1, 1, 0.3)
        # TODO: Make the colorpalette below a bit lighter
        tableau10_colors = ['006BA4', 'FF800E', 'ABABAB', '595959', '5F9ED1', 'C85200', '898989', 'A2C8EC', 'FFBC79',
                            'CFCFCF']
    else:
        fc = 'black'
        fc_i = 'white'
        fc_n = (0, 0, 0, 0.1)
        tableau10_colors = ['006BA4', 'FF800E', 'ABABAB', '595959', '5F9ED1', 'C85200', '898989', 'A2C8EC', 'FFBC79',
                            'CFCFCF']

    prop_cycle = cycler(color=['#' + s for s in tableau10_colors])
    plt.rcParams['axes.prop_cycle'] = prop_cycle

    global c
    c = prop_cycle.by_key()['color']

    plt.rcParams.update({
        'text.color': fc,
        'axes.labelcolor': fc,
        'xtick.color': fc,
        'ytick.color': fc,
        'xtick.direction': 'in',
        'ytick.direction': 'in',
        'grid.color': fc_i,
        # 'axes.axisbelow': False,
        'grid.alpha': 0.5,
        'axes.facecolor': fc_n,
        # 'axes.grid.which': 'both',
        # 'axes.grid.axis': 'both',
        'figure.facecolor': (0, 0, 0, 0),
        # 'figure.edgecolor': 'black',
        'figure.dpi': '300',
        'savefig.facecolor': (0, 0, 0, 0),
    })

    if use_kpfonts:
        plt.rcParams.update({
            'font.family': 'serif',
            'text.usetex': True,
            'text.latex.preamble': [
                r'\usepackage{amsmath}',
                r'\usepackage{amssymb}',
                r'\usepackage{siunitx}',
                r'\usepackage[notextcomp]{kpfonts}',
            ],
        })

    else:
        plt.rcParams.update({
            'font.family': 'serif',
            'font.serif': 'Times New Roman',
            'font.sans-serif': 'Times New Roman',
            'mathtext.fontset': 'cm',
        })

