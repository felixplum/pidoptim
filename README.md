# pidoptim

Demotool

# Specs:

## DataGen:
- Generiert Testdaten fuer lineares/nichtlin. Modell
- Output: t, u, y
- R/W nach .csv (von/zu Numpy)

## ModellFitter:
- Input: Alle gemessenen u, y; Modellordnungen (Verzoegerung u_lag, y_lag)
- Output: Modellfkt, welche letzte die letzten u_k.., y_(k-1).. auf y_k abbildet
- [x] Fitting of closed loop model
- [ ] Fitting of open loop model + PID params

#### Fitted closed loop model fct.
<img src="https://devfiles.syno-iq.de/s/iArmqping92Txds/preview" width="200px"/> <img src="https://devfiles.syno-iq.de/s/SY4nYNMXiFX4kme/preview" width="200px"/> <img src="https://devfiles.syno-iq.de/s/ZGdmWAgkQRjXRip/preview" width="200px"/>

## PidFitter

- Gegeben Modell (+Integrator), Referenztrajektorien, finde Kp, Kd, Ki, sodass Kosten minimal
