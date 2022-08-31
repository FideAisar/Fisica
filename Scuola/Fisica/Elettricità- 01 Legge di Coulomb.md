---
title: Elettricità - Legge di Coulomb
author: Federico Cesari
tags: [Doc]
---
## Conduzione elettrica 
### Conduzione per strofinio
Se strofino un palloncino contro della lana questo si **elettrifica**, cioè diventa capace di *attirare* o *respingere* a se altri oggetti.
>Ipotesi di Franklin: se due corpi hanno cariche opposte si attraggono, se hanno la stessa carica si respingono

Se strofino una barra metallica su di un panno di lana il panno riceve degli elettroni dalla barra che si carica positivamente mentre il panno negativamente: **effetto triboelettrico**.
! *La carica totale rimane sempre la stessa* -> [[Elettricità- 01 Legge di Coulomb#^f8569f|Conservazione della carica elettrica]]

- Si dice **isolante** elettrico un oggetto che tende a caricarsi sempre se strofinato. Questi materiali (come la gomma) sono caratterizzati da degli **elettroni molto legati** al nucleo che tendono quindi a non liberarsi e a mantenere la carica.
- Si dicono **conduttori** i materiali come i metalli che tendono a perdere la carica. Questi hanno elettroni più liberi che tendono a far perdere la carica

### Conduzione per contatto
Se un corpo carico tocca un altro conduttore, cede parte dei suoi elettroni in eccesso caricandolo negativamente. La carica elettrica $Q$ si divide equamente nei due corpi.

Il modello microscopico di Thomson ha evidenziato che la causa dell'*elettrificazione* sono gli elettroni compresi nel sistema protoni + neutroni nel nucleo atomico ed elettroni all'esterno. Difatti quando vi è un eccesso di elettroni (aumentano le cariche negative) il corpo si dice **negativo se eccede in elettroni**, **positivo se eccede in protoni**.
! *Solo gli elettroni sono mobili (un protone non si può staccare dal nucleo)*

## Carica elettrica
#### Elettroscopio
![[Electroscope.png]]
Toccando il pomello con una barra carica positivamente le due alette (spesso d'oro) si caricano entrambe positivamente e si respingono. 
>Per misurare la carica elettrica si sceglie una carica di riferimento rilevando l'ampiezza dell'angolo formato dalle due alette e in base a quello si misurano le altre cariche.

La carica elettrica è misurata in **coulomb** $C$ (dallo scienziato Charles Augustin de Coulomb) e si definisce a partire dalla carica dell'elettrone $-e = -1.6022 \cdot 10^{-19} C$. 
La carica negativa di un coulomb è quindi:
$$
\frac{-1C}{-1.6022 \cdot 10^{-19}C} =6.2414\cdot 10^{18}
\tag{1}
$$
>La carica elettrica è sempre multiplo della carica elementare $e = 1.6022 \cdot 10^{-19} C$


### Legge della conservazione della carica elettrica
^f8569f
>Indipendentemente dai fenomeni che hanno luogo, se un sistema non scambia materia con l'ambiente, la somma algebrica di tutte le cariche presenti rimarrà uguale.

## Legge di Coulomb nel vuoto
La forza elettrica di due cariche puntiformi è: 
1. Direttamente proporzionale alla forza di ciascuna carica
2. Inversamente proporzionale al quadrato della loro distanza

### Forma vettoriale
Bisogna inizialmente definire un versore $\vec{r}$ che va dal $corpo_1$ al $corpo_2$. Così possiamo esprimere la legge di Coulomb:
$$
\vec{F}=k_{0}\frac{Q_1Q_2}{r^2}\vec{r}
\tag{2}
$$
- $Q_{1,2}$ sono le due cariche
- $r$ è la distanza tra le cariche
- $k_{0}$ è la costante di proporzionalità  (nel vuoto = $8.988 \cdot 10^{9}$ $N\cdot m^{2}/C^{2}$). Spesso è conveniente scriverla come: 

$$
\frac{1}{4\pi \varepsilon_0}
$$
con $\varepsilon_{0} = 8.854 \cdot 10^{-12} \frac{C^2}{Nm^2}$ 


 ### Forma scalare
 Il modulo della forza ha la stessa equazione della forma vettoriale **senza moltiplicarlo per il versore**. 
$$
F=k_{0}\frac{|Q_1||Q_2|}{r^2}
\tag{3}
$$
>**Principio di sovrapposizione:** La forza totale che agisce su una carica elettrica è la somma vettoriale di tutte le forze delle cariche vicine come se ognuna agisse singolarmente sulla prima carica.

## Legge di Coulomb nella materia
Ovviamente se due cariche si trovano in un mezzo materiale e non nel vuoto la forza rilevata sarà minore o maggiore in relazione se il mezzo è un isolante o un conduttore.
Tramite degli esperimenti sappiamo che **il rapporto tra la forza rilevata tra due cariche nel vuoto e la forza rilevata in un mezzo è costante**. Questa costante viene chiamata **costante dielettrica relativa**.
! *Nel vuoto la costante dielettrica relativa è uguale a* $1$

La **forza è inversamente proporzionale alla costante dielettrica** del mezzo. Essendo quest'ultima posta al denominatore, al suo aumentare la forza diminuisce.
$$
\frac{F_1}{F_{2}}= \frac{\varepsilon_2}{\varepsilon_1} 
\tag{4}
$$
$$
\varepsilon_r = \frac{F}{F_m} 
\tag{5}
$$
4. Forze in due mezzi materiali diversi
5. Forze in un mezzo materiale e nel vuoto ($\varepsilon_0 = 1$)Da

Ricaviamo quindi dalla $(5)$ l'equazione della Forza in un mezzo:

$$
F_m = \frac{F}{\varepsilon_r}
\tag{6}
$$
$$
F_m = \frac{k_0}{\varepsilon_r}\frac{|Q_1||Q_2|}{r^2} 
\tag{7}
$$
Con:
- $F$ il modulo della forza nel vuoto
- $\varepsilon_r$ costante dielettrica relativa

## Attinenze con la Forza Gravitazionale
Sia la forza elettromagnetica sia quella gravitazionale:
- Agiscono a distanza e a contatto.
- Sono **direttamente proporzionali** a una grandezza caratteristica del corpo considerato (la carica e la massa).
- Sono **inversamente proporzionali** al quadrato della distanza.
Nonostante queste somiglianze però la forza elettromagnetica, a livello microscopico, è molto più forte di quella gravitazionale. 
Per dimostrarlo poniamo in esame la forza elettrica $F_e$ e gravitazionale $F_g$ tra un elettrone "$e$" e un protone "$p$" posti a una distanza $r$.

Forza elettrica:
$$
F_e = k_{0} \frac{e^2}{r^2}
\tag{8}
$$
Forza gravitazionale: 
$$
F_{g}=G \frac{m_pm_e}{r^2}
\tag{9}
$$
Con:
- $G = 6.67 \times 10^{-11} \frac{Nm^{2}}{kg^2}$ 
- $m_p=1.67\times 10^{-27} kg$
-  $m_e=9.11\times 10^{-31} kg$

Il rapporto tra le due forze è:
$$
\frac{F_e}{F_g} = \frac{k_{0}\: e^{2}}{Gm_{p}m_{e}}=2.27 \times 10^{39}
\tag{10}
$$
La forza elettrica **supera di 39 ordini di grandezza quella gravitazionale** che diventa del tutto trascurabile.

## Esperimento di Coulomb
![[bilanciaDiTorsioneCoulomb.png]]
Si hanno due cariche: una appesa a un supporto isolante e una retta da un braccio isolante e bilanciata con un peso facendo si che possa ruotare attorno al centro. La rotazione è innescata dalla forza repulsiva $F$ presente tra le due cariche. Ruotando, la seconda carica, torce un filo appeso in cima al macchinario perpendicolare al braccio.
>Il modulo $M$ del momento elastico generato dalla torsione è direttamente proporzionale all'angolo $\theta$ formato dal filo.

$$
M = c\theta
\tag{11}
$$
In cui $c$ è un valore dato dal tipo di filo utilizzato.
Il momento però può anche essere espresso come la Forza $F$ per il braccio $b$ arrivando quindi a ricavare la forza elettrica $F$ in funzione dell'angolo $\theta$ e del braccio $b$:
$$
M = Fb = c\theta
$$
$$
Fb = c\theta
$$
$$
F = c \frac{\theta}{b} 
\tag{12}
$$
