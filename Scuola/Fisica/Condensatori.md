---
title: Condensatori
author: Federico Cesari
tags: [Doc]
date: 22.09.2022
---
# Condensatori
## Carica di un condensatore
Per caricare un conduttore bisogna compiere *lavoro* per vincere la forza repulsiva delle cariche, lo stesso vale nei condensatori. Quando una forza esterna carica un condensatore esercitando un *lavoro*, l'energia rimane immagazzinata fino a quando il condensatore non è posto in un circuito.

Utilizziamo come esempio un circuito con un generatore la cui resistenza è trascurabile:

```tikz
\usepackage{pgfplots}

\begin{document}

\begin{tikzpicture}
\begin{axis}[
    title=Carica di un condensatore,
    xlabel=$\textnormal{intervallo di tempo } t$,
    ylabel=$\textnormal{intensità di corrente }i$,
    xticklabels={0,0}, yticklabels={0,0},
    xmajorgrids=true, ymajorgrids=true,
    grid style=dashed,
    clip=false
];

\addplot[samples=40,
         domain=0:10]{exp(-x)}
    node[pos = 1]{$5RC$};

\end{axis}
\end{tikzpicture}

\end{document}
```


![400](circuito_1.svg)

La *forza elettromotrice* $f_{em}$ è la differenza di potenziale massima ai capi di un generatore elettrico sconnesso dal circuito. Questa sarà sempre maggiore della $d.d.p.$ utile a causa della resistenza interna al generatore.
![400](graficorc.svg)

