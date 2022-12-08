# **FINANCIAL MARKET SIMULATOR**

## Description

This market simulator utilizes real world stock market data from the S&P500. The two main tools divide the project into two categories:
- User Portfolio - The user is allowed to input the funds that are to be invested in the market. A few predefined strategies at different risk levels are presented. The user is asked to choose their strategy, according to their risk profile. The different stocks composing each category are available to be outputted. After, the user is asked to introduce the assets to be used. According to the different asset allocation strategies and the parameters introduced, different coefficients will distribute how the funds are to be utilized.

- Market Simulator - With the user's portfolio built, the program utilizes _Geometric Brownian Motion_. It is a continuous-time stochastic process that models the stock price of the selected stocks, based on existing historical data from the S&P500. With the user's input of the investment time, the program cycles obtaining the last stock price, showing how each stock has evolved. Return on each stock and portfolio return on investment can be consulted. The program can be executed as many times as wanted to test the different strategies with different assets.

## Installation
For the program to work, the user must download all the documents presented.

External requirements: 
- C++11
- gcc compiler version > 10.3.1

## Authors
- Aníbal Guerrero Hernández
- Ignacio Zúñiga Martinez

M.Sc. Aerospace TUM Students WS22/23

## Support
Under any technical issues or further consultation, please refer to us through:

- anibal.guerrero@tum.de
- i.zuniga@tum.de

## Roadmap
Current release: **Sprint 1**
Future releases:
- Sprint 2: The objective is to create more classes for wallet identification and improve those existing for shares and strategies. The goal is to create further interactions between the market simulation and the wallets simulation, meaning that a massive shares purchase in the market leads to a change in the value of the shares themselves. The statistical data of the market fluctuations and quotes will be registered, and may be utilized in the single investment strategies.
- Sprint 3: Code performance will be evaluated and improved. Bulk tests will be executed, finding potential bottlenecks and investigating optimization options for the code, including loop transformations, STL containers, vectorization options. Complexity analysis to find potential improvements too.

## License
As an open source project, please consult the _GNU GPLv3_ license in LICENSE.md.
