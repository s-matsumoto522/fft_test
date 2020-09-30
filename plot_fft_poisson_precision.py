import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import cm
from ipywidgets import interact, FloatSlider
import numpy as np
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
#スタイルの設定
sns.set(font="TakaoPGothic", palette="colorblind", style="whitegrid")
#datファイルからデータを読み取る
axisx, axisy = np.loadtxt("./chk_fft_poisson_precision.dat", comments='!', unpack=True)
#x, yの値の設定
x = axisx
y = axisy
log_x = [np.log(i) for i in x]
log_y = [np.log(i) for i in y]
#回帰直線
a, b = np.polyfit(log_x, log_y, 1)
print('a is :', a)
print('b is :', b)
# 1. Figureのインスタンスを生成
fig = plt.figure()
# 2. Axesのインスタンスを生成
ax = fig.add_subplot(111)
# 3. データを渡してプロット
line1, = ax.plot(x, y, label = 'poisson_precision', marker='.')
# 4. グラフタイトル, ラベル付け等
ax.set_xlim(2.0*np.pi/256, 2.0*np.pi/32)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("dX")
ax.text(0.1, 0.00022, 'slope : 2.001')
ax.set_title("check the poisson precision")
ax.legend()
# 5. グラフを描画
plt.grid(which="both")
plt.show()