import numpy as np
import numpy.linalg as LA
from inputGF import *
import sys
from multiprocessing import Process, Queue
from surfGF import *

def worker(g, klist, queue, cnt):
    surf_spectral = []
    surf_modified_spectral = []
    for kpt in klist:
        tmp, tmp1 = per_k(g,kpt)
        surf_spectral += [tmp]
        surf_modified_spectral += [tmp1]
    surf_spectral = np.array(surf_spectral)
    surf_modified_spectral = np.array(surf_modified_spectral)
    queue.put((cnt, surf_spectral, surf_modified_spectral))


if __name__=='__main__':
    process_num = int(sys.argv[1])
    g = Input()
    klist = np.array(g.klist)
    pkpts = len(klist) // process_num
    klist_list = []
    for cnt in range(process_num):
        if cnt == process_num - 1:
            klist_temp = klist[cnt * pkpts:, :]
        else:
            klist_temp = klist[cnt * pkpts:(cnt + 1) * pkpts, :]
        klist_list.append(klist_temp)
    # set up queue
    queue = Queue()
    # spawn processes
    jobs = []
    for cnt in range(process_num):
        job = Process(target=worker, args=(g, klist_list[cnt], queue, cnt))
        jobs.append(job)
        job.start()
    # get results
    results = []
    for cnt in range(process_num):
        result = queue.get()
        results.append(result)
    # join all processes
    for job in jobs:
        job.join()
    # sort results list
    results = sorted(results, key=lambda result: result[0])
    # concatenate results
    surf_spectral = np.concatenate([result[1] for result in results], axis=0)
    surf_modified_spectral = np.concatenate([result[2] for result in results], axis=0)
    # write result
    print(len(g.klist), len(g.eng_list), len(surf_spectral), len(surf_spectral[0]))
    file = open('surface.dat', 'w')
    print(g.dict['mode'], file=file)
    if g.dict['mode'] == 'band':
        for i in range(len(g.klist)):
            kpt = g.klist[i]
            for j in range(len(g.eng_list)):
                print('{0:4d}    {1:12f}    {2:10f}    {3:10f}    {4:10f}    {5:12f}    {6:12f}    {7:12f}'.format(i, g.kd[i],
                kpt[0], kpt[1], kpt[2], g.eng_list[j], surf_spectral[i][j], surf_modified_spectral[i][j]), file=file)
    elif g.dict['mode'] == 'fermi':
        for i in range(len(g.klist)):
            kpt = g.point_scale(g.klist[i],g.b)
            for j in range(len(g.eng_list)):
                print('{0:4d}    {1:12f}    {2:10f}    {3:10f}    {4:12f}    {5:12f}    {6:12f}'.format(i,
                kpt[0], kpt[1], kpt[2], g.eng_list[j], surf_spectral[i][j], surf_modified_spectral[i][j]), file=file)
    file.close()