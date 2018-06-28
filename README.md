# pidoptim

Demotool

# Specs:

## DataGen:
----(Stefan 25.06.2018)----
- Generiert Testdaten fuer lineares/nichtlin. Modell
- Output: t, u, y
- R/W nach .csv (von/zu Numpy)

## ModellFitter:
----(Felix 26.06.2018)----
- generateModel(coeff, n_order) --> generiert Modellfunktion y_k = f(y_k_-1 , .. u_k, params), welche fuer Fitting und Validierung verwendet wird
	-	Zuerst: ARX-Modell; spaeter neurales Netz o.ae.
- generateObjective(u_in, y_out, params)
 	- Wird an Solver gegeben; ruft intern modelfct. auf
- Modell als Casadiexpressions (Vorschlag Felix, to be discussed)

![alt text](https://devfiles.syno-iq.de/s/iArmqping92Txds/preview "model1")


## PidFitter

- Gegeben Modell (+Integrator), Referenztrajektorien, finde Kp, Kd, Ki, sodass Kosten minimal
