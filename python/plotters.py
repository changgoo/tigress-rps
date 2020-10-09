from startup import *
from mpl_toolkits.axes_grid1.axes_grid import ImageGrid,AxesGrid
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.ticker as ticker
from matplotlib.colors import ListedColormap

def texteffect(fontsize=12,foreground='w',linewidth=3):
    try:
        from matplotlib.patheffects import withStroke
        myeffect = withStroke(foreground=foreground,linewidth=linewidth)
        kwargs = dict(path_effects=[myeffect], fontsize=fontsize)
    except ImportError:
        kwargs = dict(fontsize=fontsize)
    return kwargs

jhcolors = [
   0.000000e+00 , 0.000000e+00 , 5.000000e-01 ,
   0.000000e+00 , 0.000000e+00 , 5.156863e-01 ,
   0.000000e+00 , 0.000000e+00 , 5.313725e-01 ,
   0.000000e+00 , 0.000000e+00 , 5.470588e-01 ,
   0.000000e+00 , 0.000000e+00 , 5.627451e-01 ,
   0.000000e+00 , 0.000000e+00 , 5.784314e-01 ,
   0.000000e+00 , 0.000000e+00 , 5.941176e-01 ,
   0.000000e+00 , 0.000000e+00 , 6.098039e-01 ,
   0.000000e+00 , 0.000000e+00 , 6.254902e-01 ,
   0.000000e+00 , 0.000000e+00 , 6.411765e-01 ,
   0.000000e+00 , 0.000000e+00 , 6.568627e-01 ,
   0.000000e+00 , 0.000000e+00 , 6.725490e-01 ,
   0.000000e+00 , 0.000000e+00 , 6.882353e-01 ,
   0.000000e+00 , 0.000000e+00 , 7.039216e-01 ,
   0.000000e+00 , 0.000000e+00 , 7.196078e-01 ,
   0.000000e+00 , 0.000000e+00 , 7.352941e-01 ,
   0.000000e+00 , 0.000000e+00 , 7.509804e-01 ,
   0.000000e+00 , 0.000000e+00 , 7.666667e-01 ,
   0.000000e+00 , 0.000000e+00 , 7.823529e-01 ,
   0.000000e+00 , 0.000000e+00 , 7.980392e-01 ,
   0.000000e+00 , 0.000000e+00 , 8.137255e-01 ,
   0.000000e+00 , 0.000000e+00 , 8.294118e-01 ,
   0.000000e+00 , 0.000000e+00 , 8.450980e-01 ,
   0.000000e+00 , 0.000000e+00 , 8.607843e-01 ,
   0.000000e+00 , 0.000000e+00 , 8.764706e-01 ,
   0.000000e+00 , 0.000000e+00 , 8.921569e-01 ,
   0.000000e+00 , 0.000000e+00 , 9.078431e-01 ,
   0.000000e+00 , 0.000000e+00 , 9.235294e-01 ,
   0.000000e+00 , 0.000000e+00 , 9.392157e-01 ,
   0.000000e+00 , 0.000000e+00 , 9.549020e-01 ,
   0.000000e+00 , 0.000000e+00 , 9.705882e-01 ,
   0.000000e+00 , 0.000000e+00 , 9.862745e-01 ,
   0.000000e+00 , 1.960784e-03 , 9.960938e-01 ,
   0.000000e+00 , 1.764706e-02 , 9.960938e-01 ,
   0.000000e+00 , 3.333333e-02 , 9.960938e-01 ,
   0.000000e+00 , 4.901961e-02 , 9.960938e-01 ,
   0.000000e+00 , 6.470588e-02 , 9.960938e-01 ,
   0.000000e+00 , 8.039216e-02 , 9.960938e-01 ,
   0.000000e+00 , 9.607843e-02 , 9.960938e-01 ,
   0.000000e+00 , 1.117647e-01 , 9.960938e-01 ,
   0.000000e+00 , 1.274510e-01 , 9.960938e-01 ,
   0.000000e+00 , 1.431373e-01 , 9.960938e-01 ,
   0.000000e+00 , 1.588235e-01 , 9.960938e-01 ,
   0.000000e+00 , 1.745098e-01 , 9.960938e-01 ,
   0.000000e+00 , 1.901961e-01 , 9.960938e-01 ,
   0.000000e+00 , 2.058824e-01 , 9.960938e-01 ,
   0.000000e+00 , 2.215686e-01 , 9.960938e-01 ,
   0.000000e+00 , 2.372549e-01 , 9.960938e-01 ,
   0.000000e+00 , 2.529412e-01 , 9.960938e-01 ,
   0.000000e+00 , 2.686275e-01 , 9.960938e-01 ,
   0.000000e+00 , 2.843137e-01 , 9.960938e-01 ,
   0.000000e+00 , 3.000000e-01 , 9.960938e-01 ,
   0.000000e+00 , 3.156863e-01 , 9.960938e-01 ,
   0.000000e+00 , 3.313725e-01 , 9.960938e-01 ,
   0.000000e+00 , 3.470588e-01 , 9.960938e-01 ,
   0.000000e+00 , 3.627451e-01 , 9.960938e-01 ,
   0.000000e+00 , 3.784314e-01 , 9.960938e-01 ,
   0.000000e+00 , 3.941176e-01 , 9.960938e-01 ,
   0.000000e+00 , 4.098039e-01 , 9.960938e-01 ,
   0.000000e+00 , 4.254902e-01 , 9.960938e-01 ,
   0.000000e+00 , 4.411765e-01 , 9.960938e-01 ,
   0.000000e+00 , 4.568627e-01 , 9.960938e-01 ,
   0.000000e+00 , 4.725490e-01 , 9.960938e-01 ,
   0.000000e+00 , 4.882353e-01 , 9.960938e-01 ,
   0.000000e+00 , 5.039216e-01 , 9.960938e-01 ,
   0.000000e+00 , 5.196078e-01 , 9.960938e-01 ,
   0.000000e+00 , 5.352941e-01 , 9.960938e-01 ,
   0.000000e+00 , 5.509804e-01 , 9.960938e-01 ,
   0.000000e+00 , 5.666667e-01 , 9.960938e-01 ,
   0.000000e+00 , 5.823529e-01 , 9.960938e-01 ,
   0.000000e+00 , 5.980392e-01 , 9.960938e-01 ,
   0.000000e+00 , 6.137255e-01 , 9.960938e-01 ,
   0.000000e+00 , 6.294118e-01 , 9.960938e-01 ,
   0.000000e+00 , 6.450980e-01 , 9.960938e-01 ,
   0.000000e+00 , 6.607843e-01 , 9.960938e-01 ,
   0.000000e+00 , 6.764706e-01 , 9.960938e-01 ,
   0.000000e+00 , 6.921569e-01 , 9.960938e-01 ,
   0.000000e+00 , 7.078431e-01 , 9.960938e-01 ,
   0.000000e+00 , 7.235294e-01 , 9.960938e-01 ,
   0.000000e+00 , 7.392157e-01 , 9.960938e-01 ,
   0.000000e+00 , 7.549020e-01 , 9.960938e-01 ,
   0.000000e+00 , 7.705882e-01 , 9.960938e-01 ,
   0.000000e+00 , 7.862745e-01 , 9.960938e-01 ,
   0.000000e+00 , 8.019608e-01 , 9.960938e-01 ,
   0.000000e+00 , 8.176471e-01 , 9.960938e-01 ,
   0.000000e+00 , 8.333333e-01 , 9.960938e-01 ,
   0.000000e+00 , 8.490196e-01 , 9.960938e-01 ,
   0.000000e+00 , 8.647059e-01 , 9.960938e-01 ,
   0.000000e+00 , 8.803922e-01 , 9.960938e-01 ,
   0.000000e+00 , 8.960784e-01 , 9.960938e-01 ,
   0.000000e+00 , 9.117647e-01 , 9.960938e-01 ,
   0.000000e+00 , 9.274510e-01 , 9.960938e-01 ,
   0.000000e+00 , 9.431373e-01 , 9.960938e-01 ,
   0.000000e+00 , 9.588235e-01 , 9.960938e-01 ,
   0.000000e+00 , 9.745098e-01 , 9.960938e-01 ,
   0.000000e+00 , 9.901961e-01 , 9.960938e-01 ,
   5.882353e-03 , 9.941176e-01 , 9.941176e-01 ,
   2.156863e-02 , 9.784314e-01 , 9.784314e-01 ,
   3.725490e-02 , 9.627451e-01 , 9.627451e-01 ,
   5.294118e-02 , 9.470588e-01 , 9.470588e-01 ,
   6.862745e-02 , 9.313725e-01 , 9.313725e-01 ,
   8.431373e-02 , 9.156863e-01 , 9.156863e-01 ,
   1.000000e-01 , 9.000000e-01 , 9.000000e-01 ,
   1.156863e-01 , 8.843137e-01 , 8.843137e-01 ,
   1.313725e-01 , 8.686275e-01 , 8.686275e-01 ,
   1.470588e-01 , 8.529412e-01 , 8.529412e-01 ,
   1.627451e-01 , 8.372549e-01 , 8.372549e-01 ,
   1.784314e-01 , 8.215686e-01 , 8.215686e-01 ,
   1.941176e-01 , 8.058824e-01 , 8.058824e-01 ,
   2.098039e-01 , 7.901961e-01 , 7.901961e-01 ,
   2.254902e-01 , 7.745098e-01 , 7.745098e-01 ,
   2.411765e-01 , 7.588235e-01 , 7.588235e-01 ,
   2.568627e-01 , 7.431373e-01 , 7.431373e-01 ,
   2.725490e-01 , 7.274510e-01 , 7.274510e-01 ,
   2.882353e-01 , 7.117647e-01 , 7.117647e-01 ,
   3.039216e-01 , 6.960784e-01 , 6.960784e-01 ,
   3.196078e-01 , 6.803922e-01 , 6.803922e-01 ,
   3.352941e-01 , 6.647059e-01 , 6.647059e-01 ,
   3.509804e-01 , 6.490196e-01 , 6.490196e-01 ,
   3.666667e-01 , 6.333333e-01 , 6.333333e-01 ,
   3.823529e-01 , 6.176471e-01 , 6.176471e-01 ,
   3.980392e-01 , 6.019608e-01 , 6.019608e-01 ,
   4.137255e-01 , 5.862745e-01 , 5.862745e-01 ,
   4.294118e-01 , 5.705882e-01 , 5.705882e-01 ,
   4.450980e-01 , 5.549020e-01 , 5.549020e-01 ,
   4.607843e-01 , 5.392157e-01 , 5.392157e-01 ,
   4.764706e-01 , 5.235294e-01 , 5.235294e-01 ,
   4.921569e-01 , 5.078431e-01 , 5.078431e-01 ,
   5.078431e-01 , 5.078431e-01 , 4.921569e-01 ,
   5.235294e-01 , 5.235294e-01 , 4.764706e-01 ,
   5.392157e-01 , 5.392157e-01 , 4.607843e-01 ,
   5.549020e-01 , 5.549020e-01 , 4.450980e-01 ,
   5.705882e-01 , 5.705882e-01 , 4.294118e-01 ,
   5.862745e-01 , 5.862745e-01 , 4.137255e-01 ,
   6.019608e-01 , 6.019608e-01 , 3.980392e-01 ,
   6.176471e-01 , 6.176471e-01 , 3.823529e-01 ,
   6.333333e-01 , 6.333333e-01 , 3.666667e-01 ,
   6.490196e-01 , 6.490196e-01 , 3.509804e-01 ,
   6.647059e-01 , 6.647059e-01 , 3.352941e-01 ,
   6.803922e-01 , 6.803922e-01 , 3.196078e-01 ,
   6.960784e-01 , 6.960784e-01 , 3.039216e-01 ,
   7.117647e-01 , 7.117647e-01 , 2.882353e-01 ,
   7.274510e-01 , 7.274510e-01 , 2.725490e-01 ,
   7.431373e-01 , 7.431373e-01 , 2.568627e-01 ,
   7.588235e-01 , 7.588235e-01 , 2.411765e-01 ,
   7.745098e-01 , 7.745098e-01 , 2.254902e-01 ,
   7.901961e-01 , 7.901961e-01 , 2.098039e-01 ,
   8.058824e-01 , 8.058824e-01 , 1.941176e-01 ,
   8.215686e-01 , 8.215686e-01 , 1.784314e-01 ,
   8.372549e-01 , 8.372549e-01 , 1.627451e-01 ,
   8.529412e-01 , 8.529412e-01 , 1.470588e-01 ,
   8.686275e-01 , 8.686275e-01 , 1.313725e-01 ,
   8.843137e-01 , 8.843137e-01 , 1.156863e-01 ,
   9.000000e-01 , 9.000000e-01 , 1.000000e-01 ,
   9.156863e-01 , 9.156863e-01 , 8.431373e-02 ,
   9.313725e-01 , 9.313725e-01 , 6.862745e-02 ,
   9.470588e-01 , 9.470588e-01 , 5.294118e-02 ,
   9.627451e-01 , 9.627451e-01 , 3.725490e-02 ,
   9.784314e-01 , 9.784314e-01 , 2.156863e-02 ,
   9.941176e-01 , 9.941176e-01 , 5.882353e-03 ,
   9.960938e-01 , 9.901961e-01 , 0.000000e+00 ,
   9.960938e-01 , 9.745098e-01 , 0.000000e+00 ,
   9.960938e-01 , 9.588235e-01 , 0.000000e+00 ,
   9.960938e-01 , 9.431373e-01 , 0.000000e+00 ,
   9.960938e-01 , 9.274510e-01 , 0.000000e+00 ,
   9.960938e-01 , 9.117647e-01 , 0.000000e+00 ,
   9.960938e-01 , 8.960784e-01 , 0.000000e+00 ,
   9.960938e-01 , 8.803922e-01 , 0.000000e+00 ,
   9.960938e-01 , 8.647059e-01 , 0.000000e+00 ,
   9.960938e-01 , 8.490196e-01 , 0.000000e+00 ,
   9.960938e-01 , 8.333333e-01 , 0.000000e+00 ,
   9.960938e-01 , 8.176471e-01 , 0.000000e+00 ,
   9.960938e-01 , 8.019608e-01 , 0.000000e+00 ,
   9.960938e-01 , 7.862745e-01 , 0.000000e+00 ,
   9.960938e-01 , 7.705882e-01 , 0.000000e+00 ,
   9.960938e-01 , 7.549020e-01 , 0.000000e+00 ,
   9.960938e-01 , 7.392157e-01 , 0.000000e+00 ,
   9.960938e-01 , 7.235294e-01 , 0.000000e+00 ,
   9.960938e-01 , 7.078431e-01 , 0.000000e+00 ,
   9.960938e-01 , 6.921569e-01 , 0.000000e+00 ,
   9.960938e-01 , 6.764706e-01 , 0.000000e+00 ,
   9.960938e-01 , 6.607843e-01 , 0.000000e+00 ,
   9.960938e-01 , 6.450980e-01 , 0.000000e+00 ,
   9.960938e-01 , 6.294118e-01 , 0.000000e+00 ,
   9.960938e-01 , 6.137255e-01 , 0.000000e+00 ,
   9.960938e-01 , 5.980392e-01 , 0.000000e+00 ,
   9.960938e-01 , 5.823529e-01 , 0.000000e+00 ,
   9.960938e-01 , 5.666667e-01 , 0.000000e+00 ,
   9.960938e-01 , 5.509804e-01 , 0.000000e+00 ,
   9.960938e-01 , 5.352941e-01 , 0.000000e+00 ,
   9.960938e-01 , 5.196078e-01 , 0.000000e+00 ,
   9.960938e-01 , 5.039216e-01 , 0.000000e+00 ,
   9.960938e-01 , 4.882353e-01 , 0.000000e+00 ,
   9.960938e-01 , 4.725490e-01 , 0.000000e+00 ,
   9.960938e-01 , 4.568627e-01 , 0.000000e+00 ,
   9.960938e-01 , 4.411765e-01 , 0.000000e+00 ,
   9.960938e-01 , 4.254902e-01 , 0.000000e+00 ,
   9.960938e-01 , 4.098039e-01 , 0.000000e+00 ,
   9.960938e-01 , 3.941176e-01 , 0.000000e+00 ,
   9.960938e-01 , 3.784314e-01 , 0.000000e+00 ,
   9.960938e-01 , 3.627451e-01 , 0.000000e+00 ,
   9.960938e-01 , 3.470588e-01 , 0.000000e+00 ,
   9.960938e-01 , 3.313725e-01 , 0.000000e+00 ,
   9.960938e-01 , 3.156863e-01 , 0.000000e+00 ,
   9.960938e-01 , 3.000000e-01 , 0.000000e+00 ,
   9.960938e-01 , 2.843137e-01 , 0.000000e+00 ,
   9.960938e-01 , 2.686275e-01 , 0.000000e+00 ,
   9.960938e-01 , 2.529412e-01 , 0.000000e+00 ,
   9.960938e-01 , 2.372549e-01 , 0.000000e+00 ,
   9.960938e-01 , 2.215686e-01 , 0.000000e+00 ,
   9.960938e-01 , 2.058824e-01 , 0.000000e+00 ,
   9.960938e-01 , 1.901961e-01 , 0.000000e+00 ,
   9.960938e-01 , 1.745098e-01 , 0.000000e+00 ,
   9.960938e-01 , 1.588235e-01 , 0.000000e+00 ,
   9.960938e-01 , 1.431373e-01 , 0.000000e+00 ,
   9.960938e-01 , 1.274510e-01 , 0.000000e+00 ,
   9.960938e-01 , 1.117647e-01 , 0.000000e+00 ,
   9.960938e-01 , 9.607843e-02 , 0.000000e+00 ,
   9.960938e-01 , 8.039216e-02 , 0.000000e+00 ,
   9.960938e-01 , 6.470588e-02 , 0.000000e+00 ,
   9.960938e-01 , 4.901961e-02 , 0.000000e+00 ,
   9.960938e-01 , 3.333333e-02 , 0.000000e+00 ,
   9.960938e-01 , 1.764706e-02 , 0.000000e+00 ,
   9.960938e-01 , 1.960784e-03 , 0.000000e+00 ,
   9.862745e-01 , 0.000000e+00 , 0.000000e+00 ,
   9.705882e-01 , 0.000000e+00 , 0.000000e+00 ,
   9.549020e-01 , 0.000000e+00 , 0.000000e+00 ,
   9.392157e-01 , 0.000000e+00 , 0.000000e+00 ,
   9.235294e-01 , 0.000000e+00 , 0.000000e+00 ,
   9.078431e-01 , 0.000000e+00 , 0.000000e+00 ,
   8.921569e-01 , 0.000000e+00 , 0.000000e+00 ,
   8.764706e-01 , 0.000000e+00 , 0.000000e+00 ,
   8.607843e-01 , 0.000000e+00 , 0.000000e+00 ,
   8.450980e-01 , 0.000000e+00 , 0.000000e+00 ,
   8.294118e-01 , 0.000000e+00 , 0.000000e+00 ,
   8.137255e-01 , 0.000000e+00 , 0.000000e+00 ,
   7.980392e-01 , 0.000000e+00 , 0.000000e+00 ,
   7.823529e-01 , 0.000000e+00 , 0.000000e+00 ,
   7.666667e-01 , 0.000000e+00 , 0.000000e+00 ,
   7.509804e-01 , 0.000000e+00 , 0.000000e+00 ,
   7.352941e-01 , 0.000000e+00 , 0.000000e+00 ,
   7.196078e-01 , 0.000000e+00 , 0.000000e+00 ,
   7.039216e-01 , 0.000000e+00 , 0.000000e+00 ,
   6.882353e-01 , 0.000000e+00 , 0.000000e+00 ,
   6.725490e-01 , 0.000000e+00 , 0.000000e+00 ,
   6.568627e-01 , 0.000000e+00 , 0.000000e+00 ,
   6.411765e-01 , 0.000000e+00 , 0.000000e+00 ,
   6.254902e-01 , 0.000000e+00 , 0.000000e+00 ,
   6.098039e-01 , 0.000000e+00 , 0.000000e+00 ,
   5.941176e-01 , 0.000000e+00 , 0.000000e+00 ,
   5.784314e-01 , 0.000000e+00 , 0.000000e+00 ,
   5.627451e-01 , 0.000000e+00 , 0.000000e+00 ,
   5.470588e-01 , 0.000000e+00 , 0.000000e+00 ,
   5.313725e-01 , 0.000000e+00 , 0.000000e+00 ,
   5.156863e-01 , 0.000000e+00 , 0.000000e+00 ,
   5.000000e-01 , 0.000000e+00 , 0.000000e+00 ]

jhcmap=ListedColormap(np.array(jhcolors).reshape(256,3))
splabels = dict(ism='ISM',icm='ICM')
phlabels = dict(hot='hot',cool='cool')
Pnorm = Normalize(1,6)
Plognorm = LogNorm(10,1.e6)
Pcmap = plt.cm.plasma
norms = {'nH': Normalize(-5,1),# LogNorm(1.e-5,10),
         'temperature': Normalize(2,8),#LogNorm(100,1.e7),
         'velocity_z': SymLogNorm(100,linscale=0.5,vmin=-1000,vmax=1000,base=10),#Normalize(-500,500),
         'vB': LogNorm(vmin=1,vmax=1000),#Normalize(-500,500),
         'specific_scalar0': LogNorm(0.01,0.04),
         'Z': LogNorm(1/3,3),
         'specific_scalar3': Normalize(-2,0),
         'pok': Pnorm,
         'ram_pok_z': Pnorm,
         'mag_pok': Pnorm,
         'beta': Normalize(0,100),
         'magnetic_field_strength': Normalize(-1,2),
         'betainv': Normalize(-1,1),
         'Pram_avg': Plognorm,
         'Ptot_avg': Plognorm,
         'Pram': Plognorm,
         'Ptot': Plognorm,
         'massflux': SymLogNorm(1.e-5,vmin=-1.e-1,vmax=1.e-1,base=10),
         'momflux': Plognorm,
         'energyflux': SymLogNorm(1.e-6,vmin=-1.e-2,vmax=1.e-2,base=10),
        }
cmaps = {'nH': cma.rainforest,
         'temperature':cma.pride,
         'velocity_z': cmo.balance,
         'specific_scalar0': cma.seasons,
         'Z': cma.seasons,
         'specific_scalar3': plt.cm.cubehelix_r,
         'pok': Pcmap,
         'ram_pok_z': Pcmap,
         'mag_pok': Pcmap,
         'magnetic_field_strength': cma.horizon_r,
         'beta': cmo.delta,
         'betainv': cma.horizon_r,
         'Pram_avg': Pcmap,
         'Ptot_avg': Pcmap,
         'Pram': Pcmap,
         'Ptot': Pcmap,
         'energy': cma.seasons,
         'M3': cma.seasons,
         'vB': cma.dusk,
         'massflux': jhcmap,
         'metalflux': jhcmap,
         'momflux': Pcmap,
         'energyflux': jhcmap,
        }
variables = {
        'Pram': r'\rho v_z^2',
        'Ptot': r'P_{\rm tot}',
        'Pram_avg': r'\overline{\rho v_z^2}',
        'Ptot_avg': r'\overline{P}_{\rm tot}',
        'Z': r'Z',
        'vB': r'v_{\mathcal{B}}',
        'M3': r'\rho v_z',
        'sicm': r's_{\rm ICM}',
        'metalflux': r'\mathcal{F}_Z',
        'massflux': r'\mathcal{F}_M',
        'energyflux': r'\mathcal{F}_E',
        'momflux': r'\mathcal{F}_p',
        'metalflux_ratio': r'\mathcal{F}_Z/\mathcal{F}_{Z, {\rm in}}',
        'massflux_ratio': r'\mathcal{F}_M/\mathcal{F}_{M, {\rm in}}',
        'energyflux_ratio': r'\mathcal{F}_E/\mathcal{F}_{E, {\rm in}}',
        'momflux_ratio': r'\mathcal{F}_p/\mathcal{F}_{p, {\rm in}}',
        }
units = {
        'Ptot': r'[{\rm cm^{-3}\,K}]',
        'massflux': r'[M_\odot\,{\rm kpc^{-2}\,yr^{-1}}]',
        'momflux': r'[k_B\,{\rm cm^{-3}\,K}]',
        'energyflux': r'[{\rm 10^{51} erg\,kpc^{-2}\,yr^{-1}}]',
        }
labels = {'nH':r'$n_H$',
          'temperature':r'$T$',
          'Z':r'$Z/Z_\odot$',
          'velocity_z':r'$v_z/({\rm km/s})$',
          'specific_scalar3':r'$s_{\rm ICM}$',
          'pok':r'$P/k_B$',
          'ram_pok_z':r'$\rho v_z^2/k_B$',
          'mag_pok':r'$P_B/k_B$',
          'magnetic_field_strength':r'$B/\mu{\rm G}$',
          'beta':r'$\beta$',
          'betainv':r'$\beta^{-1}$',
         }

def plot_slices(slc,cutax='x'):
    fields = ['nH','temperature','velocity_z','Z','specific_scalar3','ram_pok_z','pok','magnetic_field_strength']
    nologs = ['specific_scalar0','Z','velocity_z','beta']
    whites = []#'nH','temperature','velocity_z','Z','specific_scalar3','magnetic_field_strength']

    nf = len(fields)
    fig=plt.figure(figsize=(2*nf,8),num=0)
    g=ImageGrid(fig,[0.05,0.1,0.85,0.75],(1,nf),label_mode='all',
                aspect=True,axes_pad=0.)

    for i,f in enumerate(fields):
        ax = g[i]
        plt.sca(ax)
        if f == 'Z':
            data = slc[cutax]['specific_scalar0']/0.02
        elif f == 'beta':
            data = slc[cutax]['pok']/slc[cutax]['mag_pok']
        elif f == 'betainv':
            data = slc[cutax]['mag_pok']/slc[cutax]['pok']
        else:
            data = slc[cutax][f]
        if not (f in nologs):
            data = np.log10(data)
#         print(f,data.min(),data.max())


        im = plt.imshow(data,extent=slc['yextent'],
                        origin='lower')
        if f in norms:
            im.set_norm(norms[f])
        if f in cmaps:
            im.set_cmap(cmaps[f])

        plt.ylim(-0.5,3.5)
        cax = inset_axes(ax,'75%','4%',loc='upper center',
                         bbox_to_anchor=[0.,0.,1.,1.1],
                         bbox_transform=ax.transAxes)

#         cax = inset_axes(ax,'75%','4%',loc='lower center',
#                          bbox_to_anchor=[0.,0.,1.,1],
#                          bbox_transform=ax.transAxes)
        cbar = plt.colorbar(cax=cax,orientation='horizontal')
        cbar.minorticks_off()
        if f == 'Z':
            cbar.set_ticks(ticker.FixedLocator([1/3,1,3]))
            cbar.ax.set_xticklabels(['1/3', '1', '3'])
        elif f == 'nH':
            cbar.set_ticks(ticker.FixedLocator([-5,-3,-1,1]))
        elif f == 'specific_scalar3':
            cbar.set_ticks(ticker.FixedLocator([-3,-2,-1,0]))
        elif f == 'temperature':
            cbar.set_ticks(ticker.FixedLocator([2,4,6,8]))
        elif 'pok' in f:
            cbar.set_ticks(ticker.FixedLocator([2,4,6]))
        elif 'velocity_z' in f:
            cbar.set_ticks([-1.e3,-1e2,1e2,1.e3])
        elif f == 'magnetic_field_strength':
            cbar.set_ticks(ticker.FixedLocator([-2,-1,0,1,2]))
        if f in labels:
            label = labels[f]
            if not (f in nologs):
                label = r'$\log($'+label+r'$)$'
            cbar.set_label(label)
#         cbar.ax.xaxis.tick_top()
        cbar.ax.xaxis.set_label_position('top')
        if (f in whites):
            set_axis_color(cbar.ax,'xy','w')
            cbar.outline.set_edgecolor('w')

    plt.setp([ax.get_yticklabels() for ax in g[1:-1]],visible=False)
    plt.setp([ax.get_xticklabels() for ax in g[1:-1]],visible=False)

    ax = g[0]
    ax.set_ylabel('z [kpc]')
    ax.set_xlabel('y [kpc]')
#     ax.xaxis.tick_top()
#     ax.xaxis.set_label_position('top')
    ax = g[-1]
    ax.set_ylabel('z [kpc]')
    ax.set_xlabel('y [kpc]')
    ax.yaxis.tick_right()
    ax.yaxis.set_label_position('right')

    for ax in g:
        ax.tick_params(axis='y',which='both',direction='inout')
        ax.tick_params(which='major',size=8)
        ax.tick_params(which='minor',size=4)
    return fig

def plot_ggsurf(prj,surf0,cbar=True,showcrit=True,cmap=cma.redshift):
    im = plt.imshow(prj['z']['data']/surf0.value,
               extent=prj['z']['bounds']/1.e3,origin='lower',
               norm=LogNorm(vmin=0.1,vmax=10),cmap=cmap)
    plt.xlabel(r'$x [{\rm kpc}]$')
    plt.ylabel(r'$y [{\rm kpc}]$')
    if cbar: plt.colorbar(im,label=r'$\Sigma_{\rm gas}/\Sigma_{\rm crit}$')
    if showcrit:
        plt.annotate(r'$\Sigma_{{\rm crit}} = {:.1f}$'.format(surf0.value) + r'$M_\odot{\rm pc^{-2}}$',
                 (0.5,1.01),xycoords='axes fraction',ha='center',va='bottom')
    return im

def get_tz_field(dc,field,sp,ph):
    zp = getattr(dc,sp).sel(phase=ph)
    if field == 'Ptot_avg':
        return zp['Ptot']/zp['A']*dc.to_pok
    elif field == 'Pth_avg':
        return zp['P']/zp['A']*dc.to_pok
    elif field == 'Pram_avg':
        return zp['Pram']/zp['A']*dc.to_pok
    elif field == 'Pmag_avg':
        return (zp['PB1'] + zp['PB2'] + zp['PB3'])/zp['A']*dc.to_pok
    elif field == 'Pimag_avg':
        return zp['Pimag']/zp['A']*dc.to_pok
    elif field == 'Ptot':
        return zp['Ptot']/dc.area*dc.to_pok
    elif field == 'Pth':
        return zp['P']/dc.area*dc.to_pok
    elif field == 'Pram':
        return zp['Pram']/dc.area*dc.to_pok
    elif field == 'Pmag':
        return (zp['PB1'] + zp['PB2'] + zp['PB3'])/dc.area*dc.to_pok
    elif field == 'Pimag':
        return zp['Pimag']/dc.area*dc.to_pok
    elif field == 'M3':
        return zp['M3']/dc.area*dc.to_flux
    elif field == 'massflux':
        return zp['massflux']/dc.area*dc.to_flux
    elif field == 'momflux':
        return zp['momflux']/dc.area*dc.to_pok
    elif field == 'energyflux':
        return zp['energyflux']/dc.area*dc.to_Eflux/1.e51
    elif field == 'Z':
        return zp['Z']/0.02
    else:
        if field in zp:
            return zp[field]
        elif field.endswith('ratio'):
            fhead = field.split('_')[0]
            return zp[fhead]/dc.picm[fhead]/dc.area

def plot_tz_all(sims,models,field,fig=None,
        boundary=False,sp='tot',ph='whole',ismcut=10,**kwargs):
    nmodel = len(models)
    if fig is None: fig=plt.figure(figsize=(10,2.5*nmodel),num=0)
    g=ImageGrid(fig,111,(nmodel,1),aspect=False,axes_pad=0.1,
                cbar_mode='single',cbar_pad=0.2,cbar_size='1%')
    iax = 0

    for m in models:
        dc = sims.dc[m]
        ax = g[iax]
        q = get_tz_field(dc,field,sp,ph)
        if (sp == 'ism') and (ismcut > 0):
            q = q.where(dc.ism['A'].sel(phase=ph) > ismcut)
        t = q.taxis*dc.to_Myr
        z = q.zaxis/1.e3
        s = dc.tot['s4'].sel(phase='whole')/dc.tot['d'].sel(phase='whole')


        plt.sca(ax)

        #if field in ['massflux','energyflux','massflux_ratio','energyflux_ratio']:
        #    _im = plt.pcolormesh(t,z,-q)
        #    if 'norm' in kwargs: _im.set_norm(kwargs['norm'])
        #    if 'cmap' in kwargs: _im.set_cmap(plt.cm.Blues)
        if sp == 'icm':
            im = plt.pcolormesh(t,z,q.where(s>1.e-15))
        else:
            im = plt.pcolormesh(t,z,q)
        if 'norm' in kwargs: im.set_norm(kwargs['norm'])
        if 'cmap' in kwargs: im.set_cmap(kwargs['cmap'])
        plt.axhline(0,ls=':',lw=3)
        plt.ylabel(r'$z [{\rm kpc}]$')

        if boundary:
            bc = (((s > 0.5)*s.zaxis).where(s>0.5)).max(dim='zaxis')/1.e3
            plt.plot(bc.taxis*dc.to_Myr,bc,color='#32cd32',lw=3,ls='-')


        plt.annotate(dc.name,(0.02,0.95),
                    xycoords='axes fraction',ha='left',va='top',
                    **texteffect('x-small'))
        iax += 1
    ax.set_xlabel(r'$t [{\rm Myr}]$')
    cbar = plt.colorbar(im,cax=g[0].cax)
    if 'label' in kwargs: cbar.set_label(kwargs['label'])
    return fig,cbar


