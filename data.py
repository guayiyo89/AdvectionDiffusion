import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#### CARGAMOS LOS DATOS DE PM
f = lambda x : (x.replace(",","."))
# PM 2.5
data = pd.read_csv("pmdatos\encinas_25.csv", sep = ",", encoding = "ISO=8859-1")
result = data.values.tolist()
pm25 = []

for a in result:
    c = -1 # default
    if pd.notna(a[2]):
        c = a[2]

    b = c
    pm25.append(b)

# PM 10
data = pd.read_csv("pmdatos\encinas_10.csv", sep = ",", encoding = "ISO=8859-1")
result = data.values.tolist()
pm10 = []

for a in result:
    c = -1 # default
    if pd.notna(a[2]):
        c = a[2]

    b = c
    pm10.append(b)

# HUMEDAD
data = pd.read_csv("datos\datos_humedad_encinas.csv", sep = ";", encoding = "ISO=8859-1", decimal = ",")
result = data.values.tolist()
datanew = np.asarray(result)

humedad = []

for a in result:
    c = -1 # Nunca se va llegar a esta T xD
    if pd.notna(a[2]):
        c = a[2]

    b = c
    humedad.append(b)

n = np.size(pm25)
m = np.size(humedad)

#### grafica PM25 / Temp
dataPM25 = []
dataPM10 = []
temp = []

for i in range(0,n):

    if pm25[i] != -1 and pm10[i] != -1 and humedad[i] != -1:
        dataPM25.append(pm25[i])
        dataPM10.append(pm10[i])
        temp.append(humedad[i])

n1 = np.size(dataPM25) #tamano muestra
m1 = np.size(temp)

dataPM25 = np.asarray(dataPM25)
dataPM10 = np.asarray(dataPM10)
temp = np.asarray(temp)

print(n1); print(m1)

def estimate_coef(x, y): 
    # number of observations/points 
    n = np.size(x) 
  
    # mean of x and y vector 
    m_x, m_y = np.mean(x), np.mean(y) 
  
    # calculating cross-deviation and deviation about x 
    SS_xy = np.sum(y*x) - n*m_y*m_x 
    SS_xx = np.sum(x*x) - n*m_x*m_x 
  
    # calculating regression coefficients 
    b_1 = SS_xy / SS_xx 
    b_0 = m_y - b_1*m_x 
  
    return(b_0, b_1)

def plot_regression_line(x, y, b): 
    # plotting the actual points as scatter plot 
    plt.scatter(x, y, color = "m", 
               marker = "o", s = 30) 
  
    # predicted response vector 
    y_pred = b[0] + b[1]*x 
  
    # plotting the regression line 
    plt.plot(x, y_pred, color = "g") 
  
    # putting labels 
    plt.xlabel('x') 
    plt.ylabel('y') 
  
    # function to show plot 
    plt.show()

b = estimate_coef(dataPM25, temp)
plot_regression_line(dataPM25, temp, b)

b2 = estimate_coef(dataPM10, temp)
plot_regression_line(dataPM10, temp, b2)