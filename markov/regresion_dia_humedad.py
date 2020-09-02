import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#### CARGAMOS LOS DATOS DE PM
f = lambda x : (x.replace(",","."))
# PM 2.5
data = pd.read_csv("pmdatos\encinaDia_25.csv", sep = ",", encoding = "ISO=8859-1")
result = data.values.tolist()
pm25 = []

for a in result:
    c = -1 # default
    if pd.notna(a[1]):
        c = a[1]

    b = c
    pm25.append(b)

# PM 10
data = pd.read_csv("pmdatos\encinaDia_10.csv", sep = ",", encoding = "ISO=8859-1")
result = data.values.tolist()
pm10 = []

for a in result:
    c = -1 # default
    if pd.notna(a[1]):
        c = a[1]

    b = c
    pm10.append(b)

# HUMEDAD
data = pd.read_csv("pmdatos\lluvia_Encina_dia.csv", encoding = "ISO=8859-1")
result = data.values.tolist()
datanew = np.asarray(result)

humedad = []

for a in result:
    c = -1 # Nunca se va llegar a esta T xD
    if pd.notna(a):
        c = a[0]

    b = c
    humedad.append(b)

n = np.size(pm25)
m = np.size(humedad)
print(m)

humedad = np.asarray(humedad)

#### grafica PM25 / Temp
dataPM25 = []
dataPM10 = []
temp = []

for i in range(0,n):

    if pm25[i] >= 0 and pm10[i] >= 0 and humedad[i] >= 0:
        dataPM25.append(pm25[i])
        dataPM10.append(pm10[i])
        temp.append(humedad[i])

n1 = np.size(dataPM25) #tamano muestra
m1 = np.size(temp)

dataPM25 = np.asarray(dataPM25)
dataPM10 = np.asarray(dataPM10)
temp = np.asarray(temp)


print(humedad)
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
    plt.plot(x, y_pred, color = "b") 
  
    # putting labels 
    plt.xlabel('lluvia') 
    plt.ylabel('PM') 
  
    # function to show plot 
    plt.show()

b = estimate_coef(temp, dataPM25)
plot_regression_line(temp, dataPM25, b)

print(b)

b1 = estimate_coef(temp, dataPM10)
plot_regression_line(temp, dataPM10, b)

print(b1)

def plot_reg2(x,y,b):
    plt.scatter(x, y, color = "m", 
               marker = "o", s = 20) 
  
    # predicted response vector 
    y_pred = b[0] + b[1]*x + b[2]*x*x
  
    # plotting the regression line 
    plt.plot(x, y_pred, color = "b") 
  
    # putting labels 
    plt.xlabel('LLUVIA') 
    plt.ylabel('PM') 
  
    # function to show plot 
    plt.show()

coeff = np.polyfit(temp, dataPM10, 2)

print(coeff)
xp = np.linspace(min(temp),max(temp),1000)
p = np.poly1d(coeff)

_ = plt.plot(temp, dataPM10, '.', xp, p(xp), '-')
plt.show()