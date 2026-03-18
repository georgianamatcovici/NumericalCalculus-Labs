import numpy as np

# metoda jacobi

def jacobi(A, eps=1e-10, kmax=100):

    n = A.shape[0]
    A = A.copy()
    U = np.eye(n)

    for k in range(kmax):

        p, q = 0, 1
        max_val = abs(A[p, q])

        for i in range(n):
            for j in range(i):
                if abs(A[i, j]) > max_val:
                    max_val = abs(A[i, j])
                    p, q = i, j

        if max_val < eps:
            break

        alpha = (A[p, p] - A[q, q]) / (2 * A[p, q])

        if alpha >= 0:
            t = -alpha + np.sqrt(alpha**2 + 1)
        else:
            t = -alpha - np.sqrt(alpha**2 + 1)

        c = 1 / np.sqrt(1 + t*t)
        s = t / np.sqrt(1 + t*t)

        for j in range(n):
            if j != p and j != q:
                apj = A[p, j]
                aqj = A[q, j]

                A[p, j] = c * apj + s * aqj
                A[q, j] = -s * apj + c * aqj

                A[j, p] = A[p, j]
                A[j, q] = A[q, j]

        app = A[p, p]
        aqq = A[q, q]
        apq = A[p, q]

        A[p, p] = app + t * apq
        A[q, q] = aqq - t * apq
        A[p, q] = 0
        A[q, p] = 0

        for i in range(n):
            uip = U[i, p]
            uiq = U[i, q]

            U[i, p] = c * uip + s * uiq
            U[i, q] = -s * uip + c * uiq

    lambdas = np.diag(A)

    return lambdas, U, A


def check_eigen(A_init, U, lambdas):

    Lambda = np.diag(lambdas)

    val = A_init @ U - U @ Lambda

    return np.linalg.norm(val)


# iteratia cholesky

def cholesky_iteration(A, eps=1e-10, kmax=100):

    A_k = A.copy()

    for k in range(kmax):

        L = np.linalg.cholesky(A_k)

        A_next = L.T @ L

        if np.linalg.norm(A_next - A_k) < eps:
            break

        A_k = A_next

    return A_k


# analiza svd (singular value decomposition)

def svd_analysis(A):

    U, S, Vt = np.linalg.svd(A)

    rank = np.sum(S > 1e-10)

    sigma_min = np.min(S)

    if sigma_min < 1e-10:
        cond_number = np.inf
    else:
        cond_number = np.max(S) / sigma_min

    return S, rank, cond_number, U, S, Vt


# pseudo-inversa svd

def pseudo_inverse_svd(A):

    U, S, Vt = np.linalg.svd(A)

    S_inv = np.zeros((Vt.shape[0], U.shape[0]))

    for i in range(len(S)):
        if S[i] > 1e-10:
            S_inv[i, i] = 1 / S[i]

    A_pinv = Vt.T @ S_inv @ U.T

    return A_pinv


# pseudo-inversa ls (least squares)

def pseudo_inverse_ls(A):

    return np.linalg.inv(A.T @ A) @ A.T


#comparatia

def compare_pseudo(A):

    AI = pseudo_inverse_svd(A)
    AJ = pseudo_inverse_ls(A)

    return np.linalg.norm(AI - AJ, 1)

def main():

    print("\nCAZUL p = n (matrice patratica):")

    A = np.array([
        [4, 1, 1],
        [1, 3, 0],
        [1, 0, 2]
    ], dtype=float)

    lambdas, U, Afinal = jacobi(A)

    print("\nMetoda Jacobi:")

    print("\nValori proprii:")
    print(lambdas)

    print("\nMatricea vectorilor proprii U:")
    print(U)

    err = check_eigen(A, U, lambdas)

    print("\nNorma ||A_init U - UΛ||:")
    print(err)

    print("\nIteratia Cholesky:")

    A_final = cholesky_iteration(A)

    print("\nMatricea finala:")
    print(A_final)

    print("\nConcluzie: Valorile de pe diagonala matricei obtinute prin iteratia Cholesky converg spre valorile proprii ale matricei initiale.\n")

    print("\nCAZUL p > n (SVD):")

    A2 = np.array([
        [1, 2],
        [3, 4],
        [5, 6]
    ], dtype=float)

    s, rank, cond, U, S, Vt = svd_analysis(A2)

    print("\nValori singulare:")
    print(s)

    print("\nRangul matricei:")
    print(rank)

    print("\nNumarul de conditie:")
    print(cond)

    AI = pseudo_inverse_svd(A2)
    AJ = pseudo_inverse_ls(A2)

    print("\nPseudo-inversa Moore-Penrose:")
    print(AI)

    print("\nPseudo-inversa Least Squares:")
    print(AJ)

    norm_val = compare_pseudo(A2)

    print("\nNorma ||AI - AJ||_1:")
    print(norm_val)


if __name__ == "__main__":
    main()
