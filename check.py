import numpy as np

def BPT(rho, n):
    """
    caculate the Bloch-Pauli Transform of a density matrix

    parameter:
        rho: density matrix

    return:
        rho_head: BPT of rho
    """
    
    X = np.array([[0,1],[1,0]])
    Y = np.array([[0,-1j],[1j,0]])
    Z = np.array([[1,0],[0,-1]])
    I = np.array([[1,0],[0,1]])

    pauli_matrices = [I,X,Y,Z]

    if n==2:
        rho_head = np.zeros(4)
        rho_head[0] = np.trace(rho @ I.conj().T).real
        rho_head[1] = np.trace(rho @ X.conj().T).real
        rho_head[2] = np.trace(rho @ Y.conj().T).real
        rho_head[3] = np.trace(rho @ Z.conj().T).real

        return rho_head
    
    elif n==4:
        rho_head = np.zeros((4,4))
        for i in range(0,4):
            for j in range(0,4):
                rho_head[i][j] = np.trace(rho @ (np.kron(pauli_matrices[i], pauli_matrices[j])).conj().T).real

        return rho_head


if __name__ == "__main__":
    n=4
    rho = np.array([[ 0  +0.j,   -0.25+0.25j,  0.  +0.j,    0.  +0.j,  ],[-0.25-0.25j,  0.5 +0.j,    0.  +0.j,    0.  +0.j,  ],[ 0.  +0.j,    0.  +0.j,    0.  +0.j,   -0.25+0.25j],[ 0.  +0.j,    0.  +0.j,   -0.25-0.25j,  0.5 +0.j  ]])
    print(BPT(rho,n))
