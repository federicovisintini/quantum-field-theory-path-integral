from scalar2d import *
import time

# Simulazione termodinamica campo scalare, energia a massa fissa
mass = 0.05
temperature = [2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 15,
               18, 20, 23, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100]
temperature = []
start = time.time()

for nt in temperature:

    # PARAMETRI
    dimension = [100, nt]
    measures = max(10 ** 6, 10 ** 5 * nt)
    termalizzazione = int(measures / 10)
    decorrel = 1

    mass2 = mass * mass
    print('temperature =', nt, '- misure =', measures / 10 ** 6, 'M')
    print('stimo servano', .5e-6 * measures * nt, 'minuti')

    simul_path = 'DATA/camposcalare/m' + str(mass) + '/nt' + str(nt) + '.dat'
    # simul_path = 'DATA/camposcalare/benchmark.dat'

    # SIMULATIONE

    # inizializzo il campo 'a freddo'
    field = np.ndarray(dimension)
    field.fill(1)

    # inizializzo il vettore che conterrà le misure (tuple di xene_mass, xene_spat, xene_temp)
    xene_mass = np.ndarray(measures)
    xene_spat = np.ndarray(measures)
    xene_temp = np.ndarray(measures)

    # termalizzazione
    for i in range(termalizzazione):
        update_heatbath(field, mass2)
        update_overrelax(field, mass2)
        update_overrelax(field, mass2)
        update_overrelax(field, mass2)
        update_overrelax(field, mass2)
        update_overrelax(field, mass2)

    # prendo 'misure' misure ogni 'decorrel' updates del reticolo
    for i in range(measures):
        for j in range(decorrel):
            update_heatbath(field, mass2)
            update_overrelax(field, mass2)
            update_overrelax(field, mass2)
            update_overrelax(field, mass2)
            update_overrelax(field, mass2)
            update_overrelax(field, mass2)

        xene_mass[i], xene_spat[i], xene_temp[i] = energy(field, mass2)

        if 100 * i % measures == 0:
            print(int(100 * i / measures + 1), '% completato')
            # e salvo la simulazione su file
            with open(simul_path, 'w') as file:
                for xene in xene_mass:
                    file.write("%f " % xene)
                file.write("\n")

                for xene in xene_spat:
                    file.write("%f " % xene)
                file.write("\n")

                for xene in xene_temp:
                    file.write("%f " % xene)

    # salvo la simulazione su file
    with open(simul_path, 'w') as file:
        for xene in xene_mass:
            file.write("%f " % xene)
        file.write("\n")

        for xene in xene_spat:
            file.write("%f " % xene)
        file.write("\n")

        for xene in xene_temp:
            file.write("%f " % xene)

print(time.time() - start)
