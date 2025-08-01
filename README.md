# clifford-density-matrix-lemma

## **Definition**
Let $\rho \in \mathcal{D}(2^n)$ represent the state of $n$ qubits. Let $\hat{\rho}$ be its BPT. We say that $\rho$ is a **Clifford density matrix** if:  

1. $\rho$ is a *scaled Hermitian projection matrix*; and  
2. All entries of $|\hat{\rho}|$ are either equal to zero or to one.

---

## **Comments**
- "Clifford density matrix" is not a standard term.  
- Let $\rho \in \mathcal{D}(2^n)$. Possibly the above Condition (2) implies the above Condition (1). If this is the case, the above definition can be simplified to listing only Condition (2). The statement in Condition (1) would then become a lemma.

---

## **About the project**
This project successfully proves the lemma in the context of single-qubit system. The result is shown in the document _single-qubit.txt_.
This project also disproves the lemma in the context of two-qubit system. The result is shown in the document _two-qubit.txt_.

--- 

## **Prerequisite**
uv installed and version >= 0.6.10


## **Usage**

First, select the number of qubit in the system in _main.py_

- Windows
```bash
.venv\Scripts\activate
uv run main.py
```
-Macos / Linux
```bash
source .venv/bin/activate
uv run main.py
```

Finally, result will be shown in _single-qubit.txt_ and _two-qubit.txt_.

