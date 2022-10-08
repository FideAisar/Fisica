---
title: Pgfplots
author: Federico Cesari
tags: [Wiki]
date: 2022/10/08
---
# Pgfplots
## 3d plots

```tikz
\usepackage{pgfplots} 
\pgfplotsset{compat=1.16} 
\begin{document} 
\begin{tikzpicture}[] 
\begin{axis}[colormap/viridis] 
\addplot3[ surf, samples=18, domain=-3:3 ] {exp(-x^2-y^2)*x}; 
\end{axis} 
\end{tikzpicture} 
\end{document} 
```

## 2d plots

```tikz
\usepackage{pgfplots} 
\pgfplotsset{compat=1.16} 
\begin{document}

\begin{tikzpicture}
	\begin{axis}[
		title={$x \exp(-x^2-y^2)$ and its gradient},
		domain=-2:2,
		view={0}{90},
		axis background/.style={fill=white},
	]
		\addplot3[contour gnuplot={number=9,
			labels=false},thick] 
				{exp(0-x^2-y^2)*x};
		\addplot3[blue,
			quiver={
			 u={exp(0-x^2-y^2)*(1-2*x^2)},
			 v={exp(0-x^2-y^2)*(-2*x*y)},
			 scale arrows=0.3,
			},
			-stealth,samples=15]
				{exp(0-x^2-y^2)*x};
	\end{axis}
\end{tikzpicture}


\end{document}
```
---

`

