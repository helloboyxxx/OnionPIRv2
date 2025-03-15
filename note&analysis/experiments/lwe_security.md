```python
params = LWE.Parameters (
    n=2048,
    q=1152921504606683137, # 60
    Xs=ND.UniformMod(3),
    Xe=ND.DiscreteGaussian(stddev=3.19), tag="ONION_2048",
)
print(params)
LWE.primal_usvp(params, red_shape_model="gsa")
```

> ```
> LWEParameters(n=2048, q=1152921504606683137, Xs=D(σ=0.82), Xe=D(σ=3.19), m=+Infinity, tag='ONION_2048')
> rop: ≈2^113.8, red: ≈2^113.8, δ: 1.004975, β: 288, d: 4056, tag: usvp
> ```



```python
params = LWE.Parameters (
    n=4096,
    q=1152921504606683137 * 1152921504606748673, # 60, 60
    Xs=ND.UniformMod (3),
    Xe=ND.DiscreteGaussian(stddev=3.19), tag="ONION_4096",
)
print(params)
LWE.primal_usvp(params, red_shape_model="gsa")
```

> ```
> LWEParameters(n=4096, q=1329227995784613643754746428306227201, Xs=D(σ=0.82), Xe=D(σ=3.19), m=+Infinity, tag='ONION_4096')
> rop: ≈2^113.7, red: ≈2^113.7, δ: 1.005021, β: 284, d: 7821, tag: usvp
> ```

