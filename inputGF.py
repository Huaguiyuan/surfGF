# Read input parameters and files
import numpy as np
import matplotlib.pyplot as plt


class Input:
    def __init__(self):
        self.a = self.read_CONTCAR()
        self.b = self.get_b(self.a)
        self.weight, self.rpt, self.hamr, self.hami = self.read_hr()  # hamr dimension: num_wann *num_wann*nrpt
        self.dict = self.read_input()
        self.klist = self.gen_klist(self.dict['k-point_div_num'], self.dict['kpath'])
        m = np.shape(self.hamr)
        self.num_wann = m[0]

        m = np.shape(self.rpt)
        self.nrpt = m[0]
        self.eng_list = self.gen_eng_list(self.dict['minimum_energy'], self.dict['maximum_energy'],
                                          self.dict['energy_div_num'])
        self.kd = self.get_k_distance(self.klist)

    def get_b(self, a):  # this funciton is to get reciprocal lattice from primitive lattice
        v = np.dot(a[0], np.cross(a[1], a[2]))
        b = []
        b.append(2 * np.pi * np.cross(a[1], a[2]) / v)
        b.append(2 * np.pi * np.cross(a[2], a[0]) / v)
        b.append(2 * np.pi * np.cross(a[0], a[1]) / v)
        b = np.array(b, dtype='float')
        return b

    def read_hr(self, file=open('wannier90_hr.dat', 'r')):
        file.readline()
        num_wann = int(file.readline().split()[0])
        nrpts = int(file.readline().split()[0])
        weight = []
        for i in range(int(np.ceil(nrpts / 15.0))):
            buffer = file.readline().split()
            weight = weight + buffer

        rpt = []
        hamr = np.zeros((num_wann, num_wann, nrpts))
        hami = np.zeros((num_wann, num_wann, nrpts))

        for i in range(nrpts):
            for j in range(num_wann):
                for k in range(num_wann):
                    buffer = file.readline().split()
                    hamr[k, j, i] = float(buffer[5])
                    hami[k, j, i] = float(buffer[6])
            rpt = rpt + [buffer[0:3]]

        weight = np.array(weight, dtype='int')
        hamr = np.array(hamr, dtype='float')
        hami = np.array(hami, dtype='float')
        rpt = np.array(rpt, dtype='int')

        file.close()

        return weight, rpt, hamr, hami

    def read_CONTCAR(self, file=open('CONTCAR', 'r')):
        file.readline()
        scale = float(file.readline())
        a1 = file.readline().split()
        a2 = file.readline().split()
        a3 = file.readline().split()
        a = [a1, a2, a3]
        a = np.array(a, dtype='float') * scale

        file.close()

        return a

    def visualize_rpt(self, rpt):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(rpt[:, 0], rpt[:, 1], rpt[:, 2])
        plt.show()
        return

    def read_input(self, file=open('INPUT', 'r')):
        import yaml
        config = file.read()
        file.close()
        dict = yaml.load(config)
        return dict

    def gen_klist(self, num, kpath):  # num per line
        kx = []
        ky = []
        kz = []
        for i in range(len(kpath) - 1):
            kx += list(np.linspace(kpath[i][0], kpath[i + 1][0], num, endpoint=False))
            ky += list(np.linspace(kpath[i][1], kpath[i + 1][1], num, endpoint=False))
            kz += list(np.linspace(kpath[i][2], kpath[i + 1][2], num, endpoint=False))
        klist = [kpath[0]]
        for i in range(len(kx)):
            klist.append([kx[i], ky[i], kz[i]])

        return klist

    def gen_eng_list(self, eng_min, eng_max, num):
        eng_list = np.linspace(eng_min, eng_max, num + 1)
        return eng_list

    def get_k_distance(self, klist):
        kd = [0]
        tmp = 0
        for i in range(len(klist) - 1):
            tmp += np.sqrt((klist[i + 1][0] - klist[i][0]) ** 2 + (klist[i + 1][1] - klist[i][1]) ** 2 + (
                klist[i + 1][2] - klist[i][2]) ** 2)
            kd.append(tmp)
        return kd