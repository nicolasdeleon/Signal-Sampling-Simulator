from tkinter import *
import csv
import matplotlib, sys
matplotlib.use('TkAgg')
import math
from scipy import *
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backends._backend_tk import NavigationToolbar2Tk
from matplotlib.figure import Figure
from matplotlib.pyplot import plot,xscale,show
import matplotlib.patches as patches
from tkinter import ttk
import numpy as np
import matplotlib.pyplot as plt
import matplotlib, sys
import matplotlib.path as mpath
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
from tkinter import messagebox
from scipy import signal
import control as ctrl

class TP1:

    def __init__(self):

        #------------------------------------Set Tk inter and Grid----------------------------------

        self.root = Tk()
        self.root.title("ASSD: Entorno de Simulacion")
        self.root.resizable(False,False)
        self.ventana_izquierda = Frame(self.root)
        self.ventana_izquierda.grid()

        self.ventana_derecha = Frame(self.root)
        self.ventana_derecha.grid(row=0,column=1)

        self.ventana_inferior = Frame(self.root)
        self.ventana_inferior.grid(row=1)

        #-----------------------------------Set Check Boxes Lists and Labels------------------------------

        self.label_1 = Label(self.ventana_derecha,text="Caracteristicas de la señal:")
        self.label_1.grid(row=0,column=2,columnspan=2, sticky=W)

        self.label_2 = Label(self.ventana_derecha,text="Etapas:")
        self.label_2.grid(row=5,column=0)

        self.label_3 = Label(self.ventana_derecha, text='Tiempo o Frecuencia:')
        self.label_3.grid(row=1, column=5, sticky=W)

        SignalList = ('Seno', '3/2 Seno','Triangular',"AM Modulada")
        self.SignalInputString = StringVar()
        self.SignalInputString.set(SignalList[0])
        InputSignalMenu = OptionMenu(self.ventana_derecha, self.SignalInputString, *SignalList)
        InputSignalMenu.grid(row=1, column=2)

        self.Frecuency = Scale(self.ventana_derecha, from_=500, to=31000, resolution=10, label='Frecuency:', orient=HORIZONTAL)
        self.Frecuency.set(1500)
        self.Frecuency.grid(row=2, column=2, padx=5, pady=5)

        self.Voltage = Scale(self.ventana_derecha, from_=0, to=10, resolution=0.1, label='Vpp:', orient=HORIZONTAL)
        self.Voltage.set(5)
        self.Voltage.grid(row=3, column=2, columnspan=1, padx=5, pady=5)

        self.check_faa = IntVar()
        self.CheckFAA = Checkbutton(self.ventana_derecha, text="FAA", variable=self.check_faa)
        self.CheckFAA.grid(row=6,column=0)

        self.check_sample_and_hold = IntVar()
        self.CheckSAH = Checkbutton(self.ventana_derecha, text="Sample and Hold", variable=self.check_sample_and_hold, command=self.mostrar_config)
        self.CheckSAH.grid(row=6,column=2,columnspan=2,padx=8,pady=8)

        self.Frecuencia_SH=Scale(self.ventana_derecha, from_=1000, to=65000, resolution=50, label='Frecuencia S&h', orient=HORIZONTAL)
        self.Frecuencia_SH.set(30000)

        self.Duty_SH=Scale(self.ventana_derecha, from_=5, to=95, label='Duty cycle S&H', orient=HORIZONTAL)
        self.Duty_SH.set(50)

        self.check_llave_analog = IntVar()
        self.CheckANKEY = Checkbutton(self.ventana_derecha, text="Llave analogica", variable=self.check_llave_analog, command=self.mostrar_config_LLA)
        self.CheckANKEY.grid(row=6,column=4,columnspan=2,padx=8,pady=8)

        self.frecuencia_LLA = Scale(self.ventana_derecha, from_=1000, to=65000,resolution=50, label='Frecuencia', orient=HORIZONTAL)
        self.frecuencia_LLA.set(3000)

        self.duty_cycle_LLA = Scale(self.ventana_derecha, from_=5, to=95, label='Duty cycle', orient=HORIZONTAL)
        self.duty_cycle_LLA.set(50)

        self.check_fr = IntVar()
        self.CheckFR = Checkbutton(self.ventana_derecha, text="FR", variable=self.check_fr)
        self.CheckFR.grid(row=6,column=6,columnspan=2,padx=8,pady=8)

        self.GraphButtom = Button(self.ventana_derecha, text="Graficar", command=self.Graficar)
        self.GraphButtom.grid(row=10,column=6,columnspan=2)

        self.check_frec_graficar = IntVar()
        self.Check_Graph_Frec = Checkbutton(self.ventana_derecha, text='Frecuencia',variable=self.check_frec_graficar)
        self.Check_Graph_Frec.grid(row=2,column=5)

        self.check_tiempo_graficar = IntVar()
        self.Check_Graph_Time = Checkbutton(self.ventana_derecha, text='Tiempo',variable=self.check_tiempo_graficar)
        self.Check_Graph_Time.grid(row=2,column=6)
         

#-----------------------------------Set Graph--------------------------------------------

        graph = Canvas(self.ventana_izquierda)
        graph.grid(row=1, columnspan=1000, padx=10, pady=10)
        f = Figure()
        self.axis = f.add_subplot(111)

        self.dataPlot = FigureCanvasTkAgg(f, master=graph)
        self.dataPlot.draw()
        self.dataPlot.get_tk_widget().pack(side=TOP, fill=BOTH, expand=True)
        self.nav = NavigationToolbar2Tk(self.dataPlot, self.ventana_inferior)
        self.nav.update()
        self.dataPlot._tkcanvas.pack(side=TOP, expand=True)

#-----------------------------RootMainLoop----------------------------------------------------

        self.root.mainloop()

#------------------------------Funcion Graficar-------
#Decidir si en Frecuencia o en Timepo. Toma valores de checkboxes de arriba para verificar esto
#Ambos seleccionados (F o T) o ambos sin seleccionar no hace nada.
#Recibe vectores datos de SetEntry
#LLama a plot tmeo plot frec respectivamente.


    def Graficar(self):
        
        t,y,NuevoN,T = self.SetEntry()
       
        if((self.check_frec_graficar.get()==0) and  (self.check_tiempo_graficar.get()==1)):
            taux = np.linspace(0.0, len(t)*T, len(t)) 
            self.PlotInTime(t,y,0)
            
        if((self.check_tiempo_graficar.get()==0) and (self.check_frec_graficar.get()==1)):
            self.PlotInFrecuency(NuevoN,T,y)  
        
        """
        Linea de codigo pensada para guardar lso vectores resultantes del grafico en un csv dentro de la carpeta del programa.
        Ojo que cada vez que se toca graficar se solapa el csv si no se cambia el nombre del .csv
        np.savetxt("Medicion17.2.csv", np.column_stack((t,y)), delimiter=",")
        """

#-------------Funcion Set Entry----------------------------------
#Toma los datos de los check boxes y de las Lists para saber que graficar y como hacerlo.
#Genera las funciones correspondientes adecuadas a la lista de funciones seleccionada
#Genera un vector de tiempo adecuado a la frecuencia de dicha señal

    def SetEntry(self):

        N = 100000 

        f = self.Frecuency.get() 
        T = 1.0 / (1000*f) #Por cada senoidal q meto me toma 1000 puntos por frecuencia, es decir muestrea 10 veces mas rapido q la señal siempre
        t = np.linspace(0.0, N*T, N*10) #toma N puntos entre (0; N/1000f)
        distanceBetweenSamples=t[1]-t[0] 

        FrecuenciaLLA = self.frecuencia_LLA.get()
        FSnH = self.Frecuencia_SH.get()

        H1 = signal.TransferFunction([((21442.0*2.0*pi)**2.0)],[1,(21442.0*2.0*pi)/6.04,(21442*2.0*pi)**2])
        H2 = signal.TransferFunction([(14150.0*2.0*pi)**2.0],[1,(14150.0*2.0*pi)/0.41,(14150.0*2.0*pi)**2])
        H3 = signal.TransferFunction([(18636.0*2.0*pi)**2.0],[1,(18636.0*2.0*pi)/1.84,(18636.0*2.0*pi)**2])
        H4 = signal.TransferFunction([(10247.0*2.0*pi)**2.0],[1,(10247.0*2.0*pi)/0.54,(10247.0*2.0*pi)**2])
        
                                                ###---ACA DEFINO TODAS MIS FUNCIONES---####
        
        if(self.SignalInputString.get() == 'Seno'):
            y = np.sin(f * 2.0*np.pi*t)
            y = y*self.Voltage.get()

        elif(self.SignalInputString.get() == 'Triangular'):
            y = signal.sawtooth(2 * np.pi * f * t,0.5)
            y = y*self.Voltage.get()


        elif(self.SignalInputString.get() == '3/2 Seno'):
            y,time = self.threeHalfsSine(t,1/f)
            y = y*self.Voltage.get()
            t=time
       
        elif(self.SignalInputString.get() == 'AM Modulada'):
            y = 0.5*np.cos(2*np.pi*t*1.8*f)+0.5*np.cos(2*np.pi*t*2.2*f)+np.cos(2*np.pi*t*2*f)
            y = y*self.Voltage.get()
            self.subNyquistFrequency(3000, 600)

        else:
            print("Error")

                                                ###---ACA DEFINO TODOS MIS MODULOS---###

        #FAA y FRR
       
        if(self.check_faa.get()):  
            t,y = self.FFA(H1,H2,H3,H4,y,t) #FFA me divide el vector 1.1            
            y = y[(int) (len(y)/1.7):] #recorte para ver de sacar mi transitorio que se me arma de que mi señal no viene desde menos inf
            t= t[(int) (len(t)/1.7):]
            t=t-t[0]     
    
        if(self.check_sample_and_hold.get()):
           y = self.sampleAndHold(FSnH,y,distanceBetweenSamples, self.Duty_SH.get()/100)
          

        if(self.check_llave_analog.get()):
           y = self.analogKey(FrecuenciaLLA,y,distanceBetweenSamples,self.duty_cycle_LLA.get()/100)
                
        if(self.check_fr.get()):
            t=t-t[0]
            t,y = self.FFA(H1,H2,H3,H4,y,t)
           

        return t,y,(int) (len(y)/10),T

    #3/2 sine
    
    def threeHalfsSine(self,time,period):
        period=period*(2/3)  #periodo senoidal
        distance=time[1]-time[0]    #distancia entre valores de la funcion seno. Seria el periodo utilizado para la señal continua
        auxTime=np.arange(time[0],4*time[len(time)-1]+time[0],distance)  #base de tiempo mayor para la senoidal de referencia
        auxSine= np.sin(auxTime*2*np.pi*(1/period))     #seno de referencia
        timeBase = time[len(time)-1]-time[0]
        amountOfPeriods = timeBase/period
        puntitosPorPeriodo=int(period/distance)  #cantidad de puntitos que entran en un periodo
        #print(puntitosPorPeriodo)
        puntitosPorPeriodoDeTMS=int(puntitosPorPeriodo*(3/2))
        #print(puntitosPorPeriodoDeTMS)
        totalRange = int(puntitosPorPeriodoDeTMS*amountOfPeriods)
        #print(totalRange)
        auxCounter=puntitosPorPeriodoDeTMS
        newFunc = np.zeros(totalRange)                #en este vector ira la amplitud de 3/2 seno, lo hago del doble de largo para que me sorbren lugares
        signChange=-1                             #contador de la cantidad de veces que hay un cambio de positivo a negativo
        state=1
        newFunc[0]=0
        for x in range(1,totalRange):
                    
            if( ( (np.sign(auxSine[x-1])!=-1) & (np.sign(auxSine[x])==-1) ) | ( (np.sign(auxSine[x-1])!=1) & (np.sign(auxSine[x])==1) ) ):         #hubo un cambio de signo
                signChange=signChange+1                                          #aumento el contador
                if(signChange%3==0):        #cuando hay dos cambios de signo de + a - tengo que invertir la funcion seno
                    state=-state                         #esta variable me indica eso
                

            if(state==1):
                newFunc[x]=-auxSine[x]
            else:
                newFunc[x]=auxSine[x]
                
        newTime=linspace(time[0], time[len(time)-1]*(3/2), len(newFunc))
        
        return newFunc, newTime

           

    #AnalogKey
    def analogKey(self,sampleFreq, signalVector, distance, dutyCycle):
        auxCounter=1/sampleFreq     #trabajo con periodo
        outputSignal = signalVector
        openKeyCondition=auxCounter*dutyCycle
        #bottomLimit=auxCounter-topLimit
        for x in range(len(signalVector)):
            if(auxCounter>=openKeyCondition):
                outputSignal[x]=0
                auxCounter=auxCounter-distance    
            else:
                if(auxCounter<=0):
                    auxCounter=1/sampleFreq         #reseteo el periodo
                outputSignal[x]=signalVector[x]
                auxCounter=auxCounter-distance

                
        return outputSignal 

    #Sample&Hold
    def sampleAndHold(self,sampleFreq, signalVector, distance, dutyCycle):
        outputSignal = signalVector
        auxCounter=1/sampleFreq
        sampleCondition=auxCounter*dutyCycle
        for x in range(len(signalVector)):
            if(auxCounter>=sampleCondition):
                holdValue=signalVector[x]
                outputSignal[x]=signalVector[x]
                auxCounter=auxCounter-distance            
            else:
                if(auxCounter<=0):
                    auxCounter=1/sampleFreq         #reseteo el periodo
                outputSignal[x]=holdValue
                auxCounter=auxCounter-distance
        return outputSignal

    def subNyquistFrequency (self,Fc,B):
        print("Rango de frecuencias subNyquist con AM modulada a 1.5k Hz")
        m=0
        Fsmax=0
        Fsmin=0
        while(Fsmax>=Fsmin):
            m=m+1
            Fsmax=(2*Fc-B)/m
            Fsmin=(2*Fc+B)/(m+1)
            print ("Para m = ", m, ": Frecuencia maxima =", Fsmax, "Frecuencia minima", Fsmin)
        m=m-1 #me quedo con el ultimo que cumple
        if(m==0):
            print("no hay frecuencia sub nyquist con la que se pueda samplear")
        Fsmax=(2*Fc-B)/m
        Fsmin=(2*Fc+B)/(m+1)
        print("Frecuencia de sampleo subNyquist para máximo m:", (Fsmax + Fsmin)/2)

        #ploteo de funcion en tiempo
    def PlotInTime(self,t,y,tipo):
        if(tipo==0):
            self.axis.clear()
            self.axis.set_title('TP1 ASSD')
            self.axis.set_aspect('auto',adjustable='box')
            self.axis.set_xlabel("Tiempo en s")
            self.axis.set_ylabel("Tensión V")
            self.axis.grid()
            self.axis.plot(t,y)
            self.dataPlot.draw()
        elif (tipo==1):
            plt.scatter(t,y)
        else:
            pass    

    def PlotInFrecuency(self,N,T,y):
        self.axis.clear()
        self.axis.set_title('TP1 ASSD')
        self.axis.set_aspect('auto',adjustable='box')
        self.axis.set_xlabel("Frecuencia en Hz")
        self.axis.set_ylabel("Tensión en dB")
        self.axis.grid()
        yfb = fft(y)
        fb = np.linspace(0.0, 1/(2.0*T),(int) (N/2)) #recorto por dos mi vector de frecuencias porque solo voy a estar con las positivas
        self.axis.loglog(fb[1:N//2], 2.0/N * np.abs(yfb[1:N//2]), '-b')
        self.dataPlot.draw()


    #Filtro FFA
    def FFA(self,H1,H2,H3,H4,y,t):
        (T1,y1,x2) = signal.lsim(H1, y, t, X0=0, interp=0)
        (T3,y2,x2) = signal.lsim(H2, y1, t, X0=0, interp=0)
        (T4,y3,x2) = signal.lsim(H3, y2, t, X0=0, interp=0)
        (T5,YAfterFilter,x2) = signal.lsim(H4, y3, t, X0=0, interp=0)
        return T5, -YAfterFilter      


    def FFR(self,H1,H2,H3,H4,y,t):
        (T2,y1,x2) = signal.lsim(H1, y, t, X0=0, interp=0)
        (T3,y2,x2) = signal.lsim(H2, y1, T2, X0=0, interp=0)
        (T4,y3,x2) = signal.lsim(H3, y2, T3, X0=0, interp=0)
        (T5,YAfterFilter,x2) = signal.lsim(H4, y3, T4, X0=0, interp=0)
        return T5, YAfterFilter        
 
    
    def mostrar_config(self):
        if (self.check_sample_and_hold.get()):
            self.Frecuencia_SH.grid(row=7, column=2, padx=8, pady=8)
            self.Duty_SH.grid(row=8,column=2,padx=8,pady=8)
        else:
            self.Frecuencia_SH.grid_remove()
            self.Duty_SH.grid_remove()

    def mostrar_config_LLA(self):
        if (self.check_llave_analog.get()):
            self.frecuencia_LLA.grid(row=7, column=4, padx=8, pady=8)
            self.duty_cycle_LLA.grid(row=8, column=4, padx=8, pady=8)
        else:
            self.frecuencia_LLA.grid_remove()
            self.duty_cycle_LLA.grid_remove()


if __name__ == "__main__":
    ex = TP1()
    ex.subNyquistFrequency(3000,600)