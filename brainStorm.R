gamma = 0.9

smallR = 1
smallT = 1
smallGa = gamma ^ smallT

largeR = 2
largeT = 2
largeGa = gamma ^ largeT

preRtimes = seq(0, 1, by = 0.2)

smallSum = smallR * (1 - smallGa ^ 10) / (1 - smallGa)
largeSum = largeR * (1 - largeGa ^ 5) / (1 - largeGa)

