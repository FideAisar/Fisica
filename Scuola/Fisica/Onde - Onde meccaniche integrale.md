---
title: Onde meccaniche integrale
author: Federico Cesari
tags: [Doc]
date: 00.04.2022
---
## Moti ondulatori
> Un'onda è una **perturbazione** che si propaga trasportando energia e quantità di moto ma **non materia**.

Un'onda si propaga a partire da un'*origine*, che nel caso di uno tsunami può trovarsi nel punto in cui è avvenuta l'eruzione vulcanica. L'onda poi (essendo meccanica) si propaga attraverso il *mezzo materiale* che è l'acqua.

### Onde trasversali e longitudinali
Se osserviamo il moto con il quale le onde si propagano attraverso il mezzo materiale riusciamo a distinguere due tipi di onde:
1. **Onde trasversali:** gli elementi del mezzo hanno un moto perpendicolare alla direzione di propagazione dell'energia. 
![[OndeTrasversali.gif.gif]]


2. **Onde longitudinali:** gli elementi del mezzo hanno un moto parallelo alla direzione di propagazione dell'energia.
![[OndeLongitudinali.gif]]
>Ogni onda che si propaga attraverso un mezzo sfruttando le sue proprietà elastiche prende il nome di **onda elastica**. *Es: corda tesa, sbarra d'acciaio, aria*

Nell'acqua possono propagarsi due tipi di onde differenti: Onde elastiche, nel caso delle onde sonore, onde non elastiche nel caso delle onde in superficie che si muovono solo grazie alla massa d'acua in movimento.

---
## Fronti d'onda
>Il **fronte d'onda** è l'insieme dei punti nei quali la variazione di grandezza al passare dell'onda ha lo stesso valore in qualunque istante. Questi oscillano concordemente, in modo tale che **la fase dell'onda passante per quei punti sia la medesima**.

Le onde però si possono propagare in tre diversi modi: 
- **In un mezzo unidimensionale**: come lungo una corda.
- **In un mezzo bidimensionale:** le onde si propagano su un piano e radialmente in tutte le direzioni. *A grande distanza dal punto d'origine le onde si possono approssimare a **rette parallele***.
- **In un mezzo tridimensionale:** le onde si propagano con la forme di tante sfere concentriche. *A grande distanza dal punto d'origine le onde si possono approssimare a **superfici piane parallele**.*
![[OndeSferiche.jpg]]
>Un'onda sferica può quindi essere descritta come un'onda piana con fronti d'onda piani e paralleli

Il **raggio d'onda** è la retta perpendicolare ai fronti d'onda.

---
## Onde periodiche
>Un'onda periodica è un'onda la cui sorgente compie un moto periodico e il cui profilo si ripete identico a distanze regolari. (Il profilo può anche essere irregolare e non sinusoidale, lo è se la sorgente oscilla in modo armonico)

### Proprietà di un'onda
1. La **lunghezza d'onda** $\lambda$ è la minima distanza dopo la quale il profilo di un'onda si ripete identico a se stesso. Rappresenta il **periodo spaziale** dell'onda.
2. L'**ampiezza** $A$ è la differenza tra il valore massimo y e il suo valore di equilibrio.
3. Il **periodo** $T$ è il tempo impiegato da y per compiere un'oscillazione completa.
4. La **frequenza** $f$ è il numero di oscillazioni compiute in $1s$.

![[LunghezzaOnda.png]]

### Velocità di propagazione
Per ogni periodo $T$ l'onda percorre la distanza di una lunghezza d'onda $\lambda$ quindi la sua velocità di propagazione ($\frac{spazio}{tempo}$) sarà:
$$
v = \frac{\lambda}{T}   
$$
>La velocità di propagazione dipende dalle proprietà del mezzo materiale.

La velocità di un'onda lungo una corda è quindi strettamente legata alla densità lineare della corda e dalla sua tensione. La formula della velocità diventa:
$$
v = \sqrt{\frac{F_T}{d_L}}
$$

---
## Onde armoniche

### Legge armonica in un punto fissato
Mentre l'onda si propaga, **un punto fissato** su di essa si muove verticalmente di un'altezza $y$ in funzione del tempo $t$.

![[SimpleHarmonicMotion.gif]]

$$
y = acos(\frac{2\pi}{T} + \upvarphi_{0})= acos(\upomega t+\upvarphi_0)
$$
- $a$ = **ampiezza** dell'onda
- $T$ = **periodo** 
- $\upomega = \frac{2\pi}{T}$ = **pulsazione** espressa in $rad/s$  
- $\frac{2\pi}{T} + \upvarphi_{0} = \upomega t+\upvarphi_0$ = argomento del coseno, è la **fase** del moto espressa in $rad$
- $\upvarphi_0$ = fase iniziale all'istante $t = 0$ 

>La **fase iniziale** è una costante che tiene conto delle condizioni iniziali dell'oscillazione.

Se consideriamo una cosinusoide di ampiezza $a$, $y$ di equilibrio = $y_0$, di funzione: $y = acos(\upomega t+\upvarphi_0)$ e poniamo $t = 0$ e $y = y_0$ allora l'equazione diventa:
$$
y = acos(\upvarphi_0)
$$
Questo dimostra che lo spostamento iniziale dipende strettamente dalla fase iniziale $\upvarphi_0$.

### Legge armonica in un istante fissato
Al variare del tempo osserviamo il profilo dell'onda, cioè il variare dell'altezza $y$ al variare dello spostamento $x$.
$$
y = acos(\frac{2\pi}{\lambda}x+\upvarphi_{0})
$$
- $\lambda$ =**lunghezza d'onda** espressa in $m$
- $x$ = **posizione x** espressa in $m$

### Funzione d'onda armonica
Dalle precedenti due funzioni che esprimono lo spostamento $y$ in funzione del tempo e dello spostamento $x$ si riesce a ricavare una formula più generale che comprende entrambi.
$$
y = acos[\frac{2\pi}{\lambda}(x-vt)+\upvarphi_0]
$$
- $v$ = velocità di propagazione dell'onda

Visto che secondo l'equazione $\frac{v}{\lambda} = \frac{1}{T}$ la formula è riscrivibile così:
$$
y = acos(\frac{2\pi}{\lambda}x-\frac{2\pi}{T}t+\upvarphi_0)
$$
Queste due equazioni possono essere utilizzate per tutti i tipi di onda, l'unico cambiamento sarebbe il significato di $y$ che in una sbarra di metallo indicherebbe la densità del metallo, per un'onda sonora il cambio di pressione, per un'onda radio il valore del campo elettrico.

---
## Interferenza
>**Principio di sovrapposizione:** quando due o più onde attraversano uno stesso mezzo producono una perturbazione uguale alla somma delle perturbazioni create singolarmente da ognuna delle onde.

Quasi tutte le onde, se non hanno ampiezze troppo elevate, soddisfano questo principio. 

### Interferenze di onde non periodiche
Gli scenari sono due: 
1. Interferenza costruttiva: quando i due effetti si sommano.
2. Interferenza distruttiva: quando i due effetti si annullano.

Se sommiamo due onde che differiscono solo per fase ($y_1=acos(\upomega t)$ e $y_2=acos(\upomega t + \upvarphi_{0})$), l'equazione risultante (prostaferesi) sarebbe:
$$
y=Acos(\upomega t+\frac{\upvarphi_0}{2})
$$
con:
$$
A= 2a\cdot cos(\frac{\upvarphi_0}{2})
$$
L'onda quindi ha una fase iniziale = $\frac{\upvarphi_0}{2}$ e l'ampiezza che dipende dal valore di $\upvarphi_0$. 
- Se $\upvarphi_{0}= 2k\pi$ con $k$ intero e quindi $cos(\frac{\upvarphi_0}{2})$ = ± 1 l'ampiezza $A$ raddoppia.
- Se $\upvarphi_{0}= (2k+1)\pi$ e $cos(\frac{\upvarphi_0}{2})$ = 0 si ha ampiezza $A=0$.

![[InterferenzaCostruttiva.gif]]

Nelle onde che si propagano su un **piano e nello spazio** invece l'interferenza è descritta per le interferenze costruttive:
$$
S_1P-S_2P = k\lambda
$$
E per le interferenze distruttive:
$$
S_1Q-S_2Q=k\lambda+\frac{1}{2}\lambda
$$
Dove:
- $S_1Q$ è la distanza dal punto di ricezione alla sorgente 1
- $S_2Q$ è la distanza dal punto di ricezione alla sorgente 2


### Sfasamento
>Lo sfasamento è dato dalla differenza $\Delta\upvarphi = \upvarphi_{2}-\upvarphi_{1}$ 

- Se lo sfasamento è nullo o è uguale a un numero intero di angoli giri. Si dice che le onde sono **in fase** = oscillano insieme nel tempo.
- Se lo sfasamento è uguale a un numero dispari di angoli piatti le onde sono in **opposizione di fase**.