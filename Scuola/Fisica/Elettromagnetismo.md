---
title: Elettromagnetismo
author: Federico Cesari
tags: [Doc]
date: 22.09.2022
---
# Elettromagnetismo
Prima di parlare di campi magnetici è necessario ricordare qualche formula e concetto riguardo il campo elettrico, l''energia potenziale e il potenziale elettrico.
## Campo elettrico
Il campo elettrico è un campo vettoriale che associa a ogni punto di uno spazio un vettore. Per misurare il campo elettrico prodotto da una o più cariche si fa uso di una *carica di prova* $q$ . La forza risultante su di essa dipenderà da ogni carica $Q$ presente nel campo e dalla sua posizione. 
La forza elettrica è definita da:
$$
|\vec{F}| = \frac{1}{4\pi \varepsilon} \frac{Qq}{r^2}  
\tag{1} \\\ \\\ \\\ \\\   

\textcolor{gray}{
\frac{Nm^2}{C^2}\frac{C^2}{m^2}
}
$$

Tuttavia, il campo elettrico esiste anche se non è presente la seconda carica $q$. Allora 
possiamo eliminarla dall'equazione definendo il modulo del campo elettrico:
$$
|\vec{E}| = \frac{|\vec{F}|}{q} \\\ \\\ \\\ \\\ \textcolor{gray}{\frac{N}{C}}
\tag{2}
$$
$$ 
|\vec{E}|=\frac{1}{4\pi \varepsilon} \frac{Q}{r^2}  \\\ \\\ \\\ \\\
\textcolor{gray}{
\frac{Nm^2}{C^2} \frac{C}{m^2} = \frac{N}{C}
}
\tag{2.1}
$$
## Energia potenziale
Sappiamo che la forza di Coulomb *è una forza conservativa*; il lavoro compiuto quindi dipende unicamente dalla posizione iniziale e da quella finale della carica.  Essendo una forza conservativa *possiamo associare ad essa un'energia potenziale*.

L'energia potenziale in un punto $A$ è uguale al lavoro che deve compiere la forza elettrica per spostare l'oggetto da $A$ fino all'infinito.
In un sistema di riferimento nel quale l'energia potenziale tra due cariche $Q_1$ e $Q_2$ è uguale a $0$ quando queste due sono infinitamente lontane tra loro, l'energia potenziale è uguale al lavoro compiuto dalla forza $F$ di $Q_1$ su $Q_2$ mentre $Q_2$ si sposta da $A$ a $B$ a distanza $r$.
$$
\Delta U = U_B - U_A = - W_{AB} 

\tag{3}
$$
La forza $F$ però, varia da punto a punto (essendo inversamente proporzionale al quadrato del raggio, allontanandosi dalla carica $Q_1$ diminuirà quadraticamente) quindi il lavoro si può trovare con l'uso di un integrale; a noi basta sapere che l'equazione che descrive l'energia potenziale è:
$$
U = \frac{1}{4\pi \varepsilon}\frac{Q_1 Q_2}{r}
\tag{4}
$$
## Potenziale elettrico
Il potenziale elettrico è strettamente legato al campo elettrico. Il campo elettrico su un oggetto $A$ infatti è una misura che dipende dall'oggetto carico e dalla distanza di questo da $A$. Queste stesse proprietà sono attribuite al *potenziale elettrico* che però è una *grandezza scalare*.
Per definire il potenziale elettrico infatti si svincola l'energia potenziale dalla sua dipendenza  con la carica di prova dividendo l'equazione per quest'ultima:
$$
V = \frac{1}{4\pi \varepsilon}\frac{Q}{r} \\\ \\\ \\\ \\\ 
\textcolor{gray}{\frac{Nm}{C} = V}
\tag{5}
$$
In un punto $A$ il potenziale sarà uguale a $V = \frac{U_A}{q}$, in $B$ sarà uguale a $V = \frac{U_B}{q}$. Definiamo quindi la *differenza di potenziale* $\Delta V$ come:
$$
\Delta V = \frac{\Delta U}{q} = -\frac{W_{AB}}{q} \\\ \\\ \\\ \\\
\textcolor{gray}{\frac{J}{C} = V} 
\tag{6}
$$
### Potenziale in un campo elettrico uniforme
Si comporta come l'energia potenziale gravitazionale in cui l'energia potenziale è uguale alla forza peso per l'altezza:
$$
U_g = mgh = Fh
$$
In questo caso però la massa è la carica $q$ e l'altezza è la distanza $z$
$$
U_e = qEz = F_e z 

$$
$$
V = \frac{U}{q} = Ez
\tag{7}
$$
---

## Carica di un condensatore
Per caricare un conduttore bisogna compiere lavoro perché le cariche già presenti su di esso oppongono resistenza via via che ne vengono aggiunte. L'equazione che descrive questo processo è:
$$
W_c =  \frac{1}{2}Q\Delta V
\tag{8}
$$
Quando il condensatore è completamente carico, la differenza di potenziale tra le due armature è uguale alla forza elettromotrice del generatore: $\Delta V =  f_{em}$ . L'equazione diventa quindi: 
$$
W_c =  \frac{1}{2}Qf_{em}
\tag{8.1}
$$
e visto che sappiamo che il lavoro compiuto dal generatore è $W_g = Qf_{em}$ (dalla $(6)$) allora *il lavoro di carica è la metà del lavoro compiuto dal generatore*. Questo ci fa capire che nel processo di carica di un condensatore metà dell'energia è utilizzata sotto forma di lavoro per portare le cariche sull'armatura, l'altra metà invece è dissipata per effetto Joule.
 
Il valore dell'intensità $i$ in funzione del tempo invece è: 
$$
i=\frac{f_{em}}{R}e^{-\frac{t}{RC}}
\tag{9}
$$
dove $RC$ è il *tempo caratteristico del circuito*.




---
Supponiamo di avere da una parte una calamita (o *dipolo magnetico*) e dall'altra due cariche unite tra loro, una positiva e una negativa. La calamita sarà indirizzata verso il polo nord terrestre in quanto il polo magnetico di questa ne è attratto. 
Pur sembrando uguali, la calamita (con i due suoi poli) e il gruppo di due cariche hanno caratteristiche diverse.  Se divido la calamita infatti, il polo attratto dal nord non sfuggirà per raggiungerlo ma diventerà una seconda calamita con due poli opposti.

Oersted si accorge che un ago magnetico sente la presenza di un filo attraverso il quale scorre una corrente elettrica. Deduce quindi
