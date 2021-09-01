from gaugeU1 import *
import time

# Simulazione topologia campo gauge, energia a massa fissa

# PARAMETRI
dimensioni = range(5, 12)
# dimensioni = [8]
measures = 1100 * 1000
decorrel = 1
termalizzazione = int(measures / 10)
delta = np.pi  # bisogna sperimentalmente trovare l'acc che da tau più piccoli
gT = 4  # rapporto fra costante di accoppiamento teoria e temperatura

# Simulation
start = time.time()
accettanza = []

for dim in dimensioni:
    # Setup
    Q = []  # lista che conterrà le misure della carica topologica
    acc = 0  # accettanza metropolis
    nx = dim
    nt = dim

    beta = nx * nt / gT ** 2  # simulazione a costante accopp. g^2 fissata
    # dichiariamo la dimensione del reticolo
    print('running nx =', nx, ' nt =', nt)
    print('ci metto =', measures * nx * nt / 600000, 'secondi')

    # Creazione delle variavili di link
    space = np.ndarray([nx, nt])
    tempo = np.ndarray([nx, nt])

    alpha = random.random() - 0.5
    for ix in range(nx):
        for it in range(nt):
            space[ix, it] = 2 * np.pi * alpha
            tempo[ix, it] = 2 * np.pi * alpha

    # start the simulation
    for i in range(termalizzazione):
        space, tempo, tmp = update_metropolis(
            space, tempo, delta, beta, nx, nt)
        # space, tempo = update_overrelax(space, tempo, nx, nt)
        # space, tempo = update_overrelax(space, tempo, nx, nt)

    for i in range(measures):
        for j in range(decorrel):
            space, tempo, tmp = update_metropolis(
                space, tempo, delta, beta, nx, nt)
            # space, tempo = update_overrelax(space, tempo, nx, nt)
            # space, tempo = update_overrelax(space, tempo, nx, nt)

        Q.append(int(round(suscett(space, tempo, nx, nt))))

        if 100 * i % measures == 0:
            print(int(100 * i / measures + 1), '% completato')

        acc += tmp

    print('accettanza metropolis = ', acc / (measures * decorrel))
    accettanza.append(acc / (measures * decorrel))

    # salvo la simulazione su file
    simul_path = 'DATA/U1gaugefield/dim' + str(dim) + '.dat'
    # simul_path = 'DATA/U1gaugefield_benchmark.dat'
    with open(simul_path, 'w') as file:
        for sushi in Q:
            file.write("%d " % sushi)

print('tempo impiegato =', time.time() - start, 'sec')
# print('accettanza:', accettanza)
