import numpy as np
import pandas as pd
import matplotlib.pyplot as plt



#### CARGAMOS LOS DATOS DE PM
f = lambda x : (x.replace(",","."))
# PM 10
data = pd.read_csv("visual/PM10_data.csv", sep = ",", encoding = "ISO=8859-1", decimal = ".")
result = data.values.tolist()
pm10 = []

cell = 24

for a in result:
    value = a[cell] * pow(10,9)
    pm10.append(value)


# PM 25
data = pd.read_csv("visual/PM25_data.csv", sep = ",", encoding = "ISO=8859-1", decimal = ".")
result = data.values.tolist()
pm25 = []

for a in result:
    value = a[cell] * pow(10,9)
    pm25.append(value)

largo = len(pm10)
x = np.arange(0, largo)

largo2 = len(pm25)
x2 = np.arange(0, largo2)

fig, axs = plt.subplots(2, sharex = True)
fig.suptitle('Concentraciones')
axs[0].plot(x, pm10)
axs[0].set_title('PM 10')
axs[0].set_ylabel('Concentracion (ug/m3)')
axs[1].plot(x2, pm25, 'tab:red')
axs[1].set_title('PM 2.5')
axs[1].set_ylabel('Concentracion (ug/m3)')
plt.xlabel('Numero de Horas')
plt.show()

#################################################################################################
#### CARGAMOS LOS DATOS DE PM
f = lambda x : (x.replace(",","."))
# PM 10
data = pd.read_csv("visual/PM10_datatest.csv", sep = ",", encoding = "ISO=8859-1", decimal = ".")
result = data.values.tolist()
pm10 = []

for a in result:
    value = a
    pm10.append(value)


# PM 25
data = pd.read_csv("visual/PM25_datatest.csv", sep = ",", encoding = "ISO=8859-1", decimal = ".")
result = data.values.tolist()
pm25 = []

# for a in result:
#     value = a
#     pm25.append(value)

# largo = len(pm10)
# x = np.arange(0, largo)
# plt.plot(x,pm10)
# plt.title("PM 10 - Conservatividad")
# plt.show()

# largo = len(pm25)
# x = np.arange(0, largo)
# plt.plot(x,pm25)
# plt.title("PM 2.5 - Conservatividad")
# plt.show()