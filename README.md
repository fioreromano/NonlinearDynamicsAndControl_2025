# 🧬 Nonlinear Dynamics and Control of a Tumor Growth Model

## 📌 Overview
This project presents an analysis and control of a nonlinear tumor growth model. The study extends a three-state differential equation system by incorporating drug absorption dynamics to enhance the realism of chemotherapy treatments.

## 🔬 Key Contributions
- **📊 Mathematical Modeling**: Formulation of a nonlinear tumor growth model, integrating the effects of chemotherapy.
- **⚖️ Equilibrium and Stability Analysis**: Identification of equilibrium points and their stability properties through linearization techniques.
- **🎯 Control Strategies**:
  - **📈 Linear Quadratic Regulator (LQR)**: Designed to drive tumor cells to extinction while preserving healthy tissue.
  - **🛡️ Sliding Mode Control (SMC)**: Implemented for robustness against model uncertainties and disturbances.
  - **🧠 Nonlinear Model Predictive Control (NMPC)**: Optimized tumor reduction while minimizing drug dosage.
- **🛠️ Robustness Analysis**: Evaluation of control strategies under parameter variations and external disturbances.
- **📊 Comparative Study**: Performance metrics comparing treatment efficiency and drug administration effectiveness.

## 📈 Results
- ✅ All controllers successfully reduced tumor size, with NMPC offering the best balance between drug minimization and treatment efficiency.
- 🛡️ SMC exhibited the highest robustness to disturbances.
- ⚙️ LQR provided a simple yet effective solution for near-equilibrium conditions.

## 🚀 Future Work
- 🏥 Extend the model to include drug resistance mechanisms.
- 🩺 Investigate patient-specific treatment personalization.

## 👩‍🎓 Author
**Fiorella Maria Romano**  
Supervised by **Prof. Mario Di Bernardo**  
📍 Scuola Politecnica e delle Scienze di Base, Ingegneria dell’Automazione e Robotica  
📅 Academic Year 2024/2025
