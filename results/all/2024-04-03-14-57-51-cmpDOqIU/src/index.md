## Pair plot 

Diagonal entries show estimates of the marginal 
densities as well as the (0.16, 0.5, 0.84) 
quantiles (dotted lines). 
Off-diagonal entries show estimates of the pairwise 
densities. 

```@raw html
<iframe src="pair_plot.png" style="height:500px;width:100%;"></iframe>
<a href="pair_plot.png"> 🔍 Full page </a>  ⏐<a href="https://sefffal.github.io/PairPlots.jl">🔗 Info </a>
```


## Trace plots 


```@raw html
<iframe src="trace_plot.png" style="height:500px;width:100%;"></iframe>
<a href="trace_plot.png"> 🔍 Full page </a>  
```


## Moments 

| **parameters** | **mean**  | **std**     | **mcse**   | **ess\_bulk** | **ess\_tail** | **rhat** | **ess\_per\_sec** |
|---------------:|----------:|------------:|-----------:|--------------:|--------------:|---------:|------------------:|
| param\_1       | 3.08537   | 0.0011397   | 7.23503e-5 | 252.15        | 504.206       | 1.00033  | missing           |
| param\_2       | 0.129818  | 0.0844345   | 0.0106568  | 63.8629       | 77.2761       | 1.01223  | missing           |
| param\_3       | -3.57424  | 0.415112    | 0.0548097  | 47.9919       | 56.2472       | 1.05697  | missing           |
| param\_4       | -0.70498  | 0.693784    | 0.0462872  | 268.295       | 86.9183       | 1.07021  | missing           |
| param\_5       | 0.812994  | 0.132876    | 0.0154557  | 53.6063       | 91.6745       | 1.04394  | missing           |
| param\_6       | -0.435046 | 1.14939     | 0.179295   | 45.6729       | 41.316        | 1.09042  | missing           |
| param\_7       | -0.994657 | 1.39955     | 0.28157    | 20.1214       | 55.9959       | 1.20204  | missing           |
| param\_8       | -1.93077  | 0.000791684 | 1.46009e-5 | 2939.53       | 3552.05       | 1.0003   | missing           |
 

```@raw html
<a href="Moments.csv">💾 CSV</a> 
```


## Cumulative traces 

For each iteration ``i``, shows the running average up to ``i``,
``\frac{1}{i} \sum_{n = 1}^{i} x_n``. 

```@raw html
<iframe src="cumulative_trace_plot.png" style="height:500px;width:100%;"></iframe>
<a href="cumulative_trace_plot.png"> 🔍 Full page </a>  
```


## Local communication barrier 

When the global communication barrier is large, many chains may 
be required to obtain tempered restarts.

The local communication barrier can be used to visualize the cause 
of a high global communication barrier. For example, if there is a 
sharp peak close to a reference constructed from the prior, it may 
be useful to switch to a [variational approximation](https://pigeons.run/dev/variational/#variational-pt).

```@raw html
<iframe src="local_barrier.png" style="height:500px;width:100%;"></iframe>
<a href="local_barrier.png"> 🔍 Full page </a>  ⏐<a href="https://pigeons.run/dev/output-pt/#Local-communication-barrier">🔗 Info </a>
```


## GCB estimation progress 

Estimate of the Global Communication Barrier (GCB) 
as a function of 
the adaptation round. 

The global communication barrier can be used 
to set the number of chains. 
The theoretical framework of [Syed et al., 2021](https://academic.oup.com/jrsssb/article/84/2/321/7056147)
yields that under simplifying assumptions, it is optimal to set the number of chains 
(the argument `n_chains` in `pigeons()`) to roughly 2Λ.

Last round estimate: ``11.339161241009801``

```@raw html
<iframe src="global_barrier_progress.png" style="height:500px;width:100%;"></iframe>
<a href="global_barrier_progress.png"> 🔍 Full page </a>  ⏐<a href="https://pigeons.run/dev/output-pt/#Global-communication-barrier">🔗 Info </a>
```


## Evidence estimation progress 

Estimate of the log normalization (computed using 
the stepping stone estimator) as a function of 
the adaptation round. 

Last round estimate: ``-27.367844403811723``

```@raw html
<iframe src="stepping_stone_progress.png" style="height:500px;width:100%;"></iframe>
<a href="stepping_stone_progress.png"> 🔍 Full page </a>  ⏐<a href="https://pigeons.run/dev/output-normalization/">🔗 Info </a>
```


## Round trips 

Number of tempered restarts  
as a function of 
the adaptation round. 

A tempered restart happens when a sample from the 
reference percolates to the target. When the reference 
supports iid sampling, tempered restarts can enable 
large jumps in the state space.

```@raw html
<iframe src="n_tempered_restarts_progress.png" style="height:500px;width:100%;"></iframe>
<a href="n_tempered_restarts_progress.png"> 🔍 Full page </a>  ⏐<a href="https://pigeons.run/dev/output-pt/#Round-trips-and-tempered-restarts">🔗 Info </a>
```


## Pigeons summary 

| **round** | **n\_scans** | **n\_tempered\_restarts** | **global\_barrier** | **global\_barrier\_variational** | **last\_round\_max\_time** | **last\_round\_max\_allocation** | **stepping\_stone** |
|----------:|-------------:|--------------------------:|--------------------:|---------------------------------:|---------------------------:|---------------------------------:|--------------------:|
| 1         | 2            | 0                         | 9.0                 | missing                          | 0.0187963                  | 640744.0                         | -1.12314e6          |
| 2         | 4            | 0                         | 4.35317             | missing                          | 0.0265403                  | 738984.0                         | -21567.1            |
| 3         | 8            | 0                         | 4.94459             | missing                          | 0.0337018                  | 428208.0                         | -14002.6            |
| 4         | 16           | 0                         | 7.97127             | missing                          | 0.0647767                  | 691232.0                         | -1001.15            |
| 5         | 32           | 0                         | 9.3428              | missing                          | 0.129562                   | 1.14685e6                        | -222.386            |
| 6         | 64           | 0                         | 10.1963             | missing                          | 0.265473                   | 2.03019e6                        | -38.541             |
| 7         | 128          | 0                         | 10.2346             | missing                          | 0.480334                   | 3.84704e6                        | -65.2633            |
| 8         | 256          | 0                         | 10.1955             | missing                          | 0.979616                   | 7.13978e6                        | -26.5775            |
| 9         | 512          | 0                         | 10.7632             | missing                          | 1.9644                     | 1.37467e7                        | -26.4203            |
| 10        | 1024         | 1                         | 11.3123             | missing                          | 3.92199                    | 2.68004e7                        | -27.7831            |
| 11        | 2048         | 11                        | 11.3607             | missing                          | 8.38295                    | 5.36352e7                        | -28.3692            |
| 12        | 4096         | 24                        | 11.3392             | missing                          | 16.1733                    | 1.04905e8                        | -27.3678            |
 

```@raw html
<a href="Pigeons_summary.csv">💾 CSV</a> ⏐<a href="https://pigeons.run/dev/output-reports/">🔗 Info </a>
```


## Pigeons inputs 

| **Keys**               | **Values**                                                                                                                                                                                                                                                            |
|-----------------------:|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| extended\_traces       | false                                                                                                                                                                                                                                                                 |
| checked\_round         | 0                                                                                                                                                                                                                                                                     |
| extractor              | nothing                                                                                                                                                                                                                                                               |
| record                 | Function[Pigeons.traces, Pigeons.round\_trip, Pigeons.log\_sum\_ratio, Pigeons.timing\_extrema, Pigeons.allocation\_extrema, Pigeons.log\_sum\_ratio, Pigeons.timing\_extrema, Pigeons.allocation\_extrema, Pigeons.round\_trip, Pigeons.energy\_ac1, Pigeons.online] |
| multithreaded          | true                                                                                                                                                                                                                                                                  |
| show\_report           | true                                                                                                                                                                                                                                                                  |
| n\_chains              | 24                                                                                                                                                                                                                                                                    |
| variational            | nothing                                                                                                                                                                                                                                                               |
| explorer               | SliceSampler(10.0, 10, 3, 1024)                                                                                                                                                                                                                                       |
| n\_chains\_variational | 0                                                                                                                                                                                                                                                                     |
| target                 | LogDensityModel for System system\_param2 of dimension 8 with fields .ℓπcallback and .∇ℓπcallback\n                                                                                                                                                                   |
| n\_rounds              | 12                                                                                                                                                                                                                                                                    |
| exec\_folder           | nothing                                                                                                                                                                                                                                                               |
| reference              | nothing                                                                                                                                                                                                                                                               |
| checkpoint             | false                                                                                                                                                                                                                                                                 |
| seed                   | 1                                                                                                                                                                                                                                                                     |
 

```@raw html
<a href="Pigeons_inputs.csv">💾 CSV</a> ⏐<a href="https://pigeons.run/dev/reference/#Pigeons.Inputs">🔗 Info </a>
```

