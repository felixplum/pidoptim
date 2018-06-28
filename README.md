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
- [x] Fitting of closed loop model
- [] Fitting of open loop model + PID params

<img src="https://devfiles.syno-iq.de/s/iArmqping92Txds/preview" alt="drawing" width="300px"/>
<img src="https://devfiles.syno-iq.de/s/SY4nYNMXiFX4kme/preview" alt="drawing" width="300px"/>
<img src="https://devfiles.syno-iq.de/s/ZGdmWAgkQRjXRip/preview" alt="drawing" width="300px"/>

## PidFitter

- Gegeben Modell (+Integrator), Referenztrajektorien, finde Kp, Kd, Ki, sodass Kosten minimal
