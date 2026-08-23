import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path



# -----------------------------------
# 参考: Physics Letters B 342(1995) 6-12, D.T. Khoa, W. von Oertzen
# a_40Ca_bdm3y1_141_7MeV
U_n_poly = [0.000949921, 0.00638314, -0.336656, 1.92118, 0.617326, 1.09765, -109.966]
U_paris_d_poly = [0.00532896, -0.143045, 1.44447, -6.64212, 13.3265, -11.0381, 3.94755, 52.221]
U_paris_e_poly = [-0.00686376, 0.187661, -1.92189, 8.92557, -18.5332, 20.6983, -7.33647, -162.125]
U_reid_d_poly = [0.00253481, -0.0669833, 0.67631, -3.19149, 6.54708, -2.6208, 2.35148, -46.1083]
U_reid_e_poly = [-0.00336784, 0.0926127, -0.958379, 4.53596, -9.6957, 10.5354, -4.88095, -62.8037]
def Extrapolater(R_mesh, R_min, R_max, f_poly_val):
    ''' 先在区间 [R_min, R_max] 实现多项式 f(为n阶多项式参数),
        然后对 R < R_min 和 R > R_max 进行外推
    '''
    in_mask = (R_mesh >= R_min) & (R_mesh <= R_max)
    left_mask = R_mesh < R_min
    right_mask = R_mesh > R_max
    f_r_mesh = np.zeros(len(R_mesh))
    f_r_mesh[in_mask] = np.polyval(f_poly_val, R_mesh[in_mask])
    # 外推
    f_r_mesh[left_mask] = f_r_mesh[in_mask][0]
    f_r_mesh[right_mask] = f_r_mesh[in_mask][-1]
    return f_r_mesh
def plot_extrapolation():
    R_mesh = np.linspace(0, 7, 1000)
    U_n = Extrapolater(R_mesh, 0.5, 8.0, U_n_poly)
    U_paris_d = Extrapolater(R_mesh, 0.5, 7.0, U_paris_d_poly)
    U_paris_e = Extrapolater(R_mesh, 0.5, 7.0, U_paris_e_poly)
    U_reid_d = Extrapolater(R_mesh, 0.5, 7.0, U_reid_d_poly)
    U_reid_e = Extrapolater(R_mesh, 0.5, 7.0, U_reid_e_poly)
    plt.plot(R_mesh, U_n, label='U_n-ref', color='blue', linestyle='-')
    plt.plot(R_mesh, U_paris_d, label='Ud-p-ref', color='red', linestyle='-')
    plt.plot(R_mesh, U_paris_e, label='Ue-p-ref', color='orange', linestyle='-')
    plt.plot(R_mesh, U_reid_d, label='Ud-r-ref', color='yellow', linestyle='-')
    plt.plot(R_mesh, U_reid_e, label='Ue-r-ref', color='green', linestyle='-')


def read_decay_file(outFilePath, colName = 'r'):
    with open(outFilePath, 'r') as f:
        lines = f.readlines()
    header = lines[0].strip().split(',')
    colnum = header.index(colName)
    colData = np.zeros(len(lines)-1)
    for i, line in enumerate(lines[1:]):
        colData[i] = float(line.strip().split(',')[colnum])
    return colData

def plot_potential(outFilePath):
    r = read_decay_file(outFilePath, colName='r')
    Ud_paris = read_decay_file(outFilePath, colName='Ud-p')
    Ue_paris = read_decay_file(outFilePath, colName='Ue-p')
    Ud_reid = read_decay_file(outFilePath, colName='Ud-r')
    Ue_reid = read_decay_file(outFilePath, colName='Ue-r')
    plt.plot(r, Ud_paris+Ue_paris, label='Un-p', linestyle='--', color='blue')
    plt.plot(r, Ud_reid+Ue_reid, label='Un-r', linestyle='--', color='blue')
    plt.plot(r, Ud_paris, label='Ud-p', linestyle='-.', color='red')
    plt.plot(r, Ue_paris, label='Ue-p', linestyle='-.', color='orange')
    plt.plot(r, Ud_reid, label='Ud-r', linestyle='-.', color='yellow')
    plt.plot(r, Ue_reid, label='Ue-r', linestyle='-.', color='green')
    # ----
    Uc_paris = read_decay_file(outFilePath, colName='Uc-p')
    Uc_reid = read_decay_file(outFilePath, colName='Uc-r')
    plt.plot(r, Uc_paris, label='Uc-p', linestyle='--', color='purple')
    plt.plot(r, Uc_reid, label='Uc-r', linestyle='-.', color='purple')

if __name__ == "__main__":
    outFilePath = Path(__file__).resolve().parents[1] / "build" / "bifold_out.csv"
    plot_extrapolation()
    plot_potential(outFilePath)
    plt.xlabel('r (fm)')
    plt.ylabel('Potential (MeV)')
    plt.title('Bifold Potentials for 40Ca at 141.7 MeV')
    plt.legend()
    plt.grid()
    plt.show()
    plt.close()
