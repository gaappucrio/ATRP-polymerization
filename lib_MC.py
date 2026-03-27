import numpy as np
import tkinter as tk
from tkinter import ttk

def Mn(DA,DB,x,L):
    MW1,MW2 = 104.14,128.17
    cont_DA = np.zeros(len(L),dtype=int)
    cont_DB = np.zeros(len(L),dtype=int)
    DA,DB = sorted(DA),sorted(DB)
    jj,kk = 0,0
    for i in range(len(L)):
        for j in range(jj,len(DA)):
            if DA[j] == x[i]:
                cont_DA[i]+=1
        jj+=cont_DA[i]
        for k in range(kk,len(DB)):
            if DB[k] == x[i]:
                cont_DB[i]+=1
        kk+=cont_DB[i]
    num = sum(x*(L+cont_DA+cont_DB))
    den = sum(L+cont_DA+cont_DB)
    Mn = num/den*(MW1+MW2)/2
    Mw = sum(x**2*(L+cont_DA+cont_DB))/num*(MW1+MW2)/2
    return Mn, Mw

def popupmsg():
    popup = tk.Tk()
    popup.wm_title("!")
    label = tk.Label(popup,text = "Simulação terminada", font = ("Verdana",12))
    label.pack(side = "top", fill = "x", pady = 10)
    B1 = tk.Button(popup, text = "Ok", command = popup.destroy)
    B1.pack()
    popup.mainloop()