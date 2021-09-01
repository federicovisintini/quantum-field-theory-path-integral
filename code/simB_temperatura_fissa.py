from scalar2d import *
import time

# Simulazione termodinamica campo scalare, energia a temperatura fissa
temperature = [25, 20, 18, 15, 12, 10, 8, 6, 5, 4, 3, 2]
start = time.time()

for nt in temperature:

    # PARAMETRI
    dimension = [80, nt]
    measures = max(1100 * 1000, 20 * 1000 * nt**2)
    termalizzazione = int(measures / 10)
    decorrel = 1
    mass = 0.1 / nt

    mass2 = mass * mass

    print('stimo servano', 0.8 * measures / 1100000 * nt / 3, 'minuti')

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

    # salvo la simulazione su file
    simul_path = 'DATA/camposcalare/T10/nt' + str(nt) + '.dat'
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
