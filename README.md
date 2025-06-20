## README

<p xmlns:cc="http://creativecommons.org/ns#" xmlns:dct="http://purl.org/dc/terms/"><span property="dct:title">Trabajo 1</span> by <span property="cc:attributionName">Aníbal Olivera</span> is licensed under <a href="https://creativecommons.org/licenses/by/4.0/?ref=chooser-v1" target="_blank" rel="license noopener noreferrer" style="display:inline-block;">CC BY 4.0<img style="height:22px!important;margin-left:3px;vertical-align:text-bottom;" src="https://mirrors.creativecommons.org/presskit/icons/cc.svg?ref=chooser-v1" alt=""><img style="height:22px!important;margin-left:3px;vertical-align:text-bottom;" src="https://mirrors.creativecommons.org/presskit/icons/by.svg?ref=chooser-v1" alt=""></a></p>

# Diffusion of Innovations with Individual Preferences: The Role of Social Reinforcement and Homophilic Ties

## About This Project

This research explores the dynamics of innovation adoption using an agent-based model that integrates individual **rational choice** with **homophilous social influence**. We aim to understand why some innovations achieve widespread success ("all-or-nothing" patterns) while others fail, by looking beyond purely structural network effects.

## Core Methodology

*   **Hybrid Agent Model:** Simulates adoption based on an innovation's intrinsic utility (`Γ`), individual preferences (`q_i`), social adoption thresholds (`τ_i`), and the scope of homophilous influence (`h`).
*   **Realistic Network Base:** Employs ATP-net (N=1000), a simulated network with socio-demographic attributes derived from the American Trends Panel.
*   **Extensive Parameter Sweep:** Analyzes over 4 million diffusion scenarios by varying:
    *   Intrinsic Innovation Utility (`Γ`)
    *   Scope of Social Influence (`h`)
    *   Mean (`μ_τ`) and Standard Deviation (`σ_τ`) of social adoption thresholds.
    *   Five distinct initial seeding strategies.

## Key Findings

Our simulations highlight that successful, widespread adoption often emerges from a critical interplay of factors:
*   The "tipping point" for mass adoption is not solely dependent on an innovation's inherent appeal (`Γ`) but is significantly modulated by the reach of social influence (`h`).
*   Increased heterogeneity in the population's social adoption thresholds (`σ_τ`) consistently promotes both higher overall adoption and the likelihood of abrupt, large-scale adoption events (phase transitions).
*   The model identifies non-structural conditions under which diffusion can be significantly blocked or rapidly accelerated.

## Repository Contents

*   `/simulation_scripts`: R scripts for the agent-based model.
*   `/analysis_scripts`: R scripts for data processing and heatmap generation.
*   `/results_data`: Raw and processed simulation output (`.rds` files).
*   `/plots`: PDF heatmaps visualizing key findings.
