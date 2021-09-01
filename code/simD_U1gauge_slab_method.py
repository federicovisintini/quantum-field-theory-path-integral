from gaugeU1 import *
import time

# Simulazione topologia campo gauge, energia a massa fissa

# PARAMETRI
dimensioni = range(22, 31)
measures = 1100 * 1000
decorrel = 1
termalizzazione = int(measures / 10)
delta = np.pi  # bisogna sperimentalmentre trovare l'acc che da tau più piccoli
gT = 4  # rapporto fra costante di accoppiamento teoria e temperatura

# Simulation
start = time.time()

for dim in dimensioni:
    # Setup
    # lista che conterrà le misure della carica topologica
    Q = np.ndarray([9, measures])
    Q_tot = []
    nx = 20
    nt = dim

    acc = 0  # accettanza metropolis
    # simulazione a costante accopp. g^2 fissata
    beta = nx * nt / gT ** 2
    # dichiariamo la dimensione del reticolo
    print('running nx =', nx, 'nt =', nt)
    print('ci metto =', measures * nx * nt / 1e6, 'secondi')

    # Creazione delle variavili di link
    space = np.ndarray([nx, nt])
    tempo = np.ndarray([nx, nt])

    # faccio tante misure indipendenti
    alpha = random.random() - 0.5
    for ix in range(nx):
        for it in range(nt):
            # space[ix, it] = 2 * np.pi * (random.random() - 1)
            # tempo[ix, it] = 2 * np.pi * (random.random() - 1)
            space[ix, it] = 2 * np.pi * alpha
            tempo[ix, it] = 2 * np.pi * alpha

    # start the simulation
    for i in range(termalizzazione):
        space, tempo, tmp = update_metropolis(
            space, tempo, delta, beta, nx, nt)
        # space, tempo = update_overrelax(space, tempo, nx, nt)

    for i in range(measures):
        for j in range(decorrel):
            space, tempo, tmp = update_metropolis(
                space, tempo, delta, beta, nx, nt)
            # space, tempo = update_overrelax(space, tempo, nx, nt)

        # misuro Q totale
        Q_tot.append(suscett(space, tempo, nx, nt))

        # e i sotto Q
        for x in range(9):
            Q[x, i] = sotto_suscett(
                space, tempo, int(nx * (x + 1) / 10), nt, nx, nt)

        if 100 * i % measures == 0:
            print(int(100 * i / measures + 1), '% completato')

        acc += tmp

    print('accettanza metropolis', acc / (measures * decorrel))

    # salvo la simulazione su file
    for x in range(9):
        simul_path = 'DATA/U1gauge_slab/dim' + \
            str(dim) + '/x' + str(x + 1) + '.dat'
        with open(simul_path, 'w') as file:
            for sushi in Q[x, :]:
                file.write("%f " % sushi)

        simul_path = 'DATA/U1gauge_slab/dim' + str(dim) + '/x0.dat'
        with open(simul_path, 'w') as file:
            for sushi in Q_tot:
                file.write("%f " % sushi)

print(time.time() - start)
