import numpy as np
import itertools
from itertools import combinations

def is_psd_matrix(mat, tol=1e-8):
    if mat.ndim != 2 or mat.shape[0] != mat.shape[1]:
        return False
    if not np.allclose(mat, mat.conj().T, atol=tol):
        return False
    eigenvalues = np.linalg.eigvalsh(mat).real
    print(eigenvalues)
    return np.all(eigenvalues >= -tol)

# Set up a BPT
def build_BPT(n,a,b,c):
    if n==2:
        rho_head = np.array([1,a,b,c])
        return rho_head
    elif n==4:
        values = [a,b,c]
        positions = [(i, j) for i in range(4) for j in range(4) if (i,j) != (0,0)]
        all_position_combinations = combinations(positions, 3)
        matrices = []
        for positions_comb in all_position_combinations:
            rho_head = np.zeros((n,n))
            rho_head[0,0] = 1
            for idx, (i, j) in enumerate(positions_comb):
                rho_head[i, j] = values[idx]
            matrices.append(rho_head)
        return matrices
    
def three_qubit_BPT(combo,f,batch_size=10000):
    print("combine")
    positions = [(i, j, k) for i in range(4) for j in range(4) for k in range(4) if (i,j,k) != (0,0,0)]
    all_position_combinations = combinations(positions, 7)
    batch = []
    print("BPT")
    for idx, positions_comb in enumerate(all_position_combinations):
        rho_head = np.zeros((4,4,4))
        rho_head[0,0,0] = 1.0
        for idx, (i, j, k) in enumerate(positions_comb):
            rho_head[i, j, k] = combo[idx]
        batch.append(rho_head)
        if len(batch) >= batch_size:
            print("overflow...")
            process_batch(batch,f)
            batch = []
    if batch:
        process_batch(batch,f)

def process_batch(batch,f):
    for rho_head in batch:
        rho = inverse_BPT(rho_head,8)
        result = is_scaled_Hermitian_projection_matrix(rho)
        statement(rho_head,result,rho,f)


def inverse_BPT(rho_head,n):
    X = np.array([[0,1],[1,0]])
    Y = np.array([[0,-1j],[1j,0]])
    Z = np.array([[1,0],[0,-1]])
    I = np.array([[1,0],[0,1]])
    pauli_matrices = [I,X,Y,Z]
    if n==2:
        rho = np.zeros((n,n), dtype=np.complex128)
        for i in range(0,4):
            rho += rho_head[i]/2*pauli_matrices[i]
    elif n==4:
        rho = np.zeros((n,n), dtype=np.complex128)
        for i in range(0,4):
            for j in range(0,4):
                rho += rho_head[i][j]*(np.kron(pauli_matrices[i]/2,pauli_matrices[j]/2))
    elif n==8:
        rho = np.zeros((n,n), dtype=np.complex128)
        for i in range(0,4):
            for j in range(0,4):
                for k in range(0,4):
                    term = rho_head[i][j][k]*(np.kron(pauli_matrices[i]/2,pauli_matrices[j]/2))
                    rho += np.kron(term,pauli_matrices[k]/2)
    return rho

def is_scaled_Hermitian_projection_matrix(rho, tol=1e-8):
    """
    check if a density matrix is a scaled Hermitian projection matrix

    parameter:
        rho: density matrix
        tol: tolerance of error

    return:
        True, c, rho
        False, None, rho
    """

    rho_square = rho @ rho

    non_zero_mask = np.abs(rho) > tol
    if np.sum(non_zero_mask) == 0:
        return False, None, rho
    
    # calculate c
    c_values = rho_square[non_zero_mask] / rho[non_zero_mask]
    c = np.mean(c_values)
    
    # check if c is nonzero，and M² ≈ c·M
    if np.abs(c) < tol:
        return False, None, rho
    if not np.allclose(rho_square, c * rho, atol=tol):
        return False, None, rho

    return True, c, rho

def norm(a,b,c):
    return np.sqrt(a*a+b*b+c*c)

def statement(rho_head,result,rho,f):
    if is_psd_matrix(rho):
        if result[0]==True:
            f.write("case:\n" + str(rho_head) + ": " + str(result[0]) + "\n")
            f.write("\n")
        elif result[0]==False:
            f.write("case:\n" + str(rho_head) + ": " + str(result[0])+"\n")
            f.write("The counterexample is:\n")
            f.write("rho_head:\n"+str(rho_head)+"\n")
            f.write("rho:\n"+str(rho)+"\n")
            f.write("\n")

def main(qubit):
    print("start")
    n = np.power(2,qubit)
    values = [-1,0,1]
    if n==2:
        with open("single-qubit.txt","w",encoding="utf-8") as f:
            for a in values:
                for b in values:
                    for c in values:
                            if norm(a,b,c) <= 1:
                                rho_head = build_BPT(n,a,b,c)
                                rho = inverse_BPT(rho_head,n)
                                result = is_scaled_Hermitian_projection_matrix(rho)
                                f.write("case "+ str(rho_head) + ": " + str(result[0]) + "\n")
    if n==4:
        with open("two-qubit.txt","w",encoding="utf-8") as f:
            all_permutations = list(itertools.product(values, repeat=3))

            unique_combinations = []
            seen = set()

            for p in all_permutations:
                sorted_p = tuple(sorted(p))
                if sorted_p not in seen:
                    seen.add(sorted_p)
                    unique_combinations.append(sorted_p)
            for combo in unique_combinations:
                a=combo[0]
                b=combo[1]
                c=combo[2]
                matrices = build_BPT(n,a,b,c)
                for rho_head in matrices:
                    rho = inverse_BPT(rho_head,n)
                    result = is_scaled_Hermitian_projection_matrix(rho)
                    statement(rho_head,result,rho,f)
    if n==8:
        with open("three-qubit.txt","w",encoding="utf-8") as f:
            print(2)
            all_permutations = list(itertools.product(values, repeat=7))

            unique_combinations = []
            seen = set()

            for p in all_permutations:
                sorted_p = tuple(sorted(p))
                if sorted_p not in seen:
                    seen.add(sorted_p)
                    unique_combinations.append(sorted_p)
            for combo in unique_combinations:
                three_qubit_BPT(combo,f)
 