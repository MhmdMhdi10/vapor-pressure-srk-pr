# vapor-pressure-srk-pr

📖 Introduction

This project calculates and plots the vapor pressure of C1, C3, C5, C7, and H₂O as a function of temperature using two thermodynamic models:

    Soave-Redlich-Kwong (SRK) Equation of State

    Peng-Robinson (PR) Equation of State

The Antoine equation provides initial saturation pressures, while fugacity corrections and the Newton-Raphson method are used to refine results.
📦 Project Structure

├── main.py/
├── Antoine_equation/
│   └── antoine.py/
├── srk_eos/
│   └── srk_eos.py/
├── peng_robinson/
│   └── peng_robinson.py/
├── newthon_raphson/
│   └── newton_raphson.py/
├── Utility/
│   ├── fv_calculator.py/
│   └── fl_calculator.py/
└── README.md/

⚙️ Installation

    Clone this repository:

git clone https://github.com/your-username/your-repo-name.git
cd your-repo-name

Install the required libraries:

    pip install numpy matplotlib scipy

🚀 How to Run

Run the main.py script:

python main.py

When prompted:

    Enter 1 to use SRK EOS.

    Enter 2 to use Peng-Robinson EOS.

The program will calculate vapor pressures and plot vapor pressure vs. temperature curves for each component.
📊 Example Output

    X-axis: Temperature (°F)

    Y-axis: Vapor Pressure (psi)

    Separate curves for each component (C1, C3, C5, C7, H₂O).

📚 References

    Antoine Equation for vapor pressure estimation.

    Soave-Redlich-Kwong (SRK) and Peng-Robinson (PR) equations of state.

    Newton-Raphson root-finding method for fugacity correction.
