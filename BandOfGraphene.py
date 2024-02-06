import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc


class BandOfGraphene:
    e2p = 0  # 2pz軌道エネルギー
    t = -3.033  # 最近接原子間のトランスファー積分
    s = 0.129  # 最近接原子間の重なり積分
    a = 2.46  # 基本格子ベクトルaの大きさ．炭素の原子間距離の√3倍
    b = (4 * np.pi) / (np.sqrt(3) * a)  # 逆格子ベクトルbの大きさ
    width = b / np.sqrt(3)  # グラフの幅

    def w(self, x: np.ndarray | float, y: np.ndarray | float) -> np.ndarray | float:
        ret = np.exp(1 * x * self.a * 1j / np.sqrt(3)) + 2 * np.exp(-1 * x * self.a * 1j / (2 * np.sqrt(3))) * np.cos(y * self.a / 2)
        return np.abs(ret)

    def E_2g(self, x: np.ndarray | float, y: np.ndarray | float, sign: str) -> np.ndarray | float:
        w = self.w(x, y)
        if sign == '+':
            return (self.e2p + self.t * w)/(1 + self.s * w)
        if sign == '-':
            return (self.e2p - self.t * w)/(1 - self.s * w)
        print('wrong sign. sign should be "+" or "-".')
        return -1

    def plot(self, n=100):
        # プロットの設定
        plt.rcParams["font.size"] = 20
        plt.rcParams["font.family"] = 'Times New Roman'
        rc('mathtext', **{'rm': 'serif', 'fontset': 'cm'})
        fig = plt.figure(figsize=(10, 6))
        ax = fig.add_subplot(111, projection='3d', facecolor='w')
        # メインの計算
        kx = np.linspace(-self.width, self.width, n)
        ky = np.linspace(-self.width, self.width, n)
        kx, ky = np.meshgrid(kx, ky)
        e2g_plus = self.E_2g(kx, ky, '+')
        e2g_minus = self.E_2g(kx, ky, '-')
        ax.plot_surface(kx, ky, e2g_plus, zorder=1, cmap='plasma', vmax=15, vmin=-10)
        ax.plot_surface(kx, ky, e2g_minus, zorder=2, cmap='plasma', vmax=15, vmin=-10)
        # 補助線，補助点
        line = np.zeros(4)
        lines_x = [self.b / 2, 0, -1 * self.b / 2, -1 * self.b / 2, 0, self.b / 2]
        lines_y = [self.width / 2, self.width, self.width / 2, -1 * self.width / 2, -1 * self.width, -1 * self.width / 2]
        labels = ['Γ', 'M', 'K']
        points_x = [0, self.b / 2, self.b / 2]
        points_y = [0, 0, self.width / 2]
        for i in range(6):
            x_tmp = np.linspace(lines_x[i - 1], lines_x[i], 4)
            y_tmp = np.linspace(lines_y[i - 1], lines_y[i], 4)
            ax.plot(x_tmp, y_tmp, line, '-', zorder=3, color='k', linewidth=2)
        for label, x, y in zip(labels, points_x, points_y):
            z = self.E_2g(x, y, '+')
            ax.plot(x, y, z, 'o', color='k', zorder=4)
            ax.text(x, y, z + 1, label, zorder=5, fontdict={'family': 'Arial'})
        # 軸ラベル
        ax.set_xlabel(r'$\rm k_{x}$')
        ax.set_ylabel(r'$\rm k_{y}$')
        ax.set_zlabel('$E$  [eV]')
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_zticks([])
        plt.show()


def main():
    b = BandOfGraphene()
    b.plot()


if __name__ == '__main__':
    main()
