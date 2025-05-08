"""
sympy helper functions
"""

import sympy  # noqa: F401
from sympy import (  # noqa: F401
    symbols,
    Matrix,
    eye,
    Function,
    diff,
    Derivative,
    cos,
    sin,
    simplify,
    MatrixBase,
    solve,
    pi,
    sqrt,
    atan,
    atan2,
    acos,
    sign,
    series,
    expand,
    integrate,
    collect,
    csc,
    tan,
    exp,
    diag,
)

# Globally useful definitions
w1, w2, w3 = symbols("omega_1,omega_2,omega_3", real=True)
t, g, m, h = symbols("t,g,m,h", real=True)
circmat = Matrix([eye(3)[2 - j, :] for j in range(3)])  # define circulant matrix
polarframe = [
    "\\hat{\\mathbf{e}}_r",
    "\\hat{\\mathbf{e}}_\\theta",
    "\\hat{\\mathbf{e}}_3",
]
polarframe_nohat = ["\\mathbf{e}_r", "\\mathbf{e}_\\theta", "\\mathbf{e}_3"]

sphericalframe = [
    "\\hat{\\mathbf{e}}_\\phi",
    "\\hat{\\mathbf{e}}_\\theta",
    "\\hat{\\mathbf{e}}_\\rho",
]
sphericalframe_nohat = ["\\mathbf{e}_\\phi", "\\mathbf{e}_\\theta", "\\mathbf{e}_\\rho"]


def gendiffvars(syms, real=True):
    """Generate symbolic variables and their derivatives

    Args:
        syms (iterable):
            List (or any iterable) of symbols to create.  Each element is either a
            string or another iterable, whose contents are:
            (variable name, [symbol name], [number of derivatives])
            If number of derivatives is not set, 2 is assumed. If symbol name is not set
            use the same name for the variable and the symbol.  Thus, an input of 'x' is
            equivalent to ('x', 'x', 2).
        real (bool):
            True for real value.

    Returns:
        tuple:
            allsyms (dict):
                All generated symbols (can be used to populate calling namespace)
            diffmap (dict):
                Differentiation map

    Examples:
        >>> # create 2nd order derivatives for theta and phi:
        >>> allsyms, diffmap = gendiffvars([('th','theta'), ('ph', 'phi')])
        >>> locals().update(allsyms)

    .. note::
        When the name of a symbol includes an underscore, the 'dot' will be placed
        preceeding the underscore such that the leading term is dotted. However, the
        corresponding variable name will always have the 'd' placed at the end,
        regardless of whether the variable name includes and underscore.  Thus, a
        variable definition like ('th_1', 'theta_1') will result in the first derivative
        being named thetadot_1 and assigned to varaible th_1d. To avoid confusion, it is
        recommended to avoid underscores in variable names - e.g., for this example to
        define the variable as ('th1', 'theta_1'), which would result in a variable th1d
        mapping to a symbol with name 'thetadot_1'.
    """

    diffmap = {}
    allsyms = {}
    for symboldef in syms:
        if isinstance(symboldef, str):
            varname = symboldef
            symname = symboldef
            nderivs = 2
        else:
            varname = symboldef[0]
            if len(symboldef) > 1:
                symname = symboldef[1]
            else:
                symname = varname
            if len(symboldef) > 2:
                nderivs = symboldef[2]
            else:
                nderivs = 2
        # generate new syms locally
        newsyms = {varname: symbols(symname, real=real)}
        for j in range(nderivs):
            # if symname includes underscore, nead to treat it differently
            if "_" in symname:
                p1, p2 = symname.split("_")
                derivname = f'{p1}{"".join(["d"] * j)}dot_{p2}'
            else:
                derivname = f'{symname}{"".join(["d"] * j)}dot'

            newsyms[f'{varname}{"".join(["d"] * (j + 1))}'] = symbols(
                derivname, real=real
            )
        locals().update(newsyms)

        # update diffmap
        for j in range(nderivs):
            diffmap[locals()[f'{varname}{"".join(["d"] * (j))}']] = locals()[
                f'{varname}{"".join(["d"] * (j + 1))}'
            ]

        # update output
        allsyms.update(newsyms)

    return allsyms, diffmap


def difftotal(expr, diffby, diffmap):
    """Take the total derivative with respect to a variable of an expression where
    dependent variables are not defined as functions of that variable

    Args:
        expr (sympy.core.*):
            Any valid sympy expression.
        diffby (sympy.core.symbol.Symbol):
            Inependent variable to differentiate with respect to
        diffmap (dict):
            Dictionary of dependent variables and their first derivatives

    Returns:
        sympy.core.*:
            Differentiated expression

    Examples:
        >>> theta, t, theta_dot = symbols("theta t theta_dot")
        >>> difftotal(cos(theta), t, {theta: theta_dot})
        -theta_dot*sin(theta)

    .. note::

        Adapted from code by Chris Wagner
        http://robotfantastic.org/total-derivatives-in-sympy.html

    """
    # Replace all symbols in the diffmap by a functional form
    fnexpr = expr.subs({s: Function(str(s))(diffby) for s in diffmap})
    # Do the differentiation
    diffexpr = diff(fnexpr, diffby)
    # Replace the Derivatives with the variables in diffmap
    derivmap = {
        Derivative(Function(str(v))(diffby), diffby): dv for v, dv in diffmap.items()
    }
    finaldiff = diffexpr.subs(derivmap)
    # Replace the functional forms with their original form
    return finaldiff.subs({Function(str(s))(diffby): s for s in diffmap})


def difftotalmat(mat, diffby, diffmap):
    """Apply total derivative element by element to a matrix

    Args:
        mat (sympy.matrices.dense.MutableDenseMatrix):
            Input matrix of expressions to differentiate
        diffby (sympy.core.symbol.Symbol):
            Inependent variable to differentiate with respect to
        diffmap (dict):
            Dictionary of dependent variables and their first derivatives

    Returns:
        sympy.matrices.dense.MutableDenseMatrix:
            Differentiated matrix (same dimensions as input matrix)

    """
    return Matrix([difftotal(x, diffby, diffmap) for x in mat]).reshape(*mat.shape)


def transportEq(vec, diffby, diffmap, omega):
    r"""Apply the transport equation to a vector. For any pair of reference frames
    :math:`\mathcal{I}` and :math:`\mathcal{B}` and any vector :math:`\mathbf{c}`:

    .. math::

        \vphantom{\frac{d}{d}}^{\mathcal{I}}\frac{\mathrm{d}}{\mathrm{d}t} \mathbf{c}=
        \vphantom{\frac{d}{d}}^{\mathcal{B}}\frac{\mathrm{d}}{\mathrm{d}t} \mathbf{c} +
        {}^\mathcal{I}\boldsymbol{\omega}^\mathcal{B} \times \mathbf{c}


    Args:
        vec (sympy.matrices.dense.MutableDenseMatrix):
            3x1 column matrix of vector components in frame :math:`\mathcal{B}`
        diffby (sympy.core.symbol.Symbol):
            Inependent variable to differentiate with respect to
        diffmap (dict):
            Dictionary of dependent variables and their first derivatives
        omega (sympy.matrices.dense.MutableDenseMatrix):
            3x1 column matrix of angular velocity components

    """

    return difftotalmat(vec, diffby, diffmap) + omega.cross(vec)


def rotMat(axis, angle):
    r"""Generate direction cosine matrix :math:`{}^\mathcal{B}C^\mathcal{A}` about
    one of the unit directions defining frame :math:`\mathcal{A}`

    Args:
        axis (int):
            1, 2, or 3 only
        angle (sympy.core.symbol.Symbol):
            angle variable

    Returns:
        sympy.matrices.dense.MutableDenseMatrix:
            The DCM
    """

    if axis == 1:
        return Matrix(
            ([1, 0, 0], [0, cos(angle), sin(angle)], [0, -sin(angle), cos(angle)])
        )
    elif axis == 2:
        return Matrix(
            ([cos(angle), 0, -sin(angle)], [0, 1, 0], [sin(angle), 0, cos(angle)])
        )
    elif axis == 3:
        return Matrix(
            ([cos(angle), sin(angle), 0], [-sin(angle), cos(angle), 0], [0, 0, 1])
        )
    else:
        return -1


def parallelAxis(I_G, r_QG, m):
    r"""Applies the parallel axis theorem to matrix of inertia I_G (where G is the
    body's center of mass) to find the matrix of inertia I_Q where the vector from G to
    Q is r_QG and the total mass of the body is m.  I_G and I_QG are assumed to be in
    components of the same (body-fixed) frame

    Args:
        I_G (sympy.matrices.dense.MutableDenseMatrix):
            Matrix of inertia about the COM in body-fixed frame components
        r_QG (sympy.matrices.dense.MutableDenseMatrix):
            3x1 matrix representation of the position of point Q with respect to G in
            components of the same body-fixed frame
        m (sympy.core.*):
            Symbol or expression of the total mass of the body

    Returns:
        sympy.matrices.dense.MutableDenseMatrix:
            resulting matrix of inertia about Q in components of the same body-fixed
            frame


    """

    return I_G + m * ((r_QG.transpose() * r_QG)[0] * eye(3) - r_QG * r_QG.transpose())


def skew(v):
    r"""Compute the skew-symmetric cross-produce equivalent matrix of a vector

    Args:
        v (sympy.matrices.dense.MutableDenseMatrix):
            3x1 matrix representation of a vector with respect to some reference frame

    Returns:
        sympy.matrices.dense.MutableDenseMatrix:
            Cross-product equivalent matrix of the vector

    """

    assert (
        hasattr(v, "__iter__") or isinstance(v, Matrix) or isinstance(v, MatrixBase)
    ) and len(v) == 3, "v must be an iterable of length 3."

    return Matrix([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])


def calcDCM(n, th):
    r"""Compute a direction cosine matrix via the Euler-Rodrigues equation.
    Evaluates the DCM :math:`{}^\mathcal{A}C^\mathcal{B}` for frames
    :math:`\mathcal{I}` and :math:`\mathcal{B}` for axis of rotation :math:`\mathbf{n}`
    and angle :math:`\theta`

    Args:
        n (sympy.matrices.dense.MutableDenseMatrix):
            3x1 matrix representation of the unit vector of the axis of rotation
        th (sympy.core.*):
            Symbol or expression for the angle of rotation

    Returns:
        sympy.matrices.dense.MutableDenseMatrix:
            The direction cosine matrix

    """

    return eye(3) * cos(th) + skew(n) * sin(th) + (1 - cos(th)) * n * n.T


def DCM2angVel(dcm, diffmap, diffby=t):
    r"""Given a direction cosine matrix :math:`{}^\mathcal{B}C^\mathcal{A}` compute
    the angular velocity :math:`{}^\mathcal{I}\boldsymbol{\omega}^\mathcal{B}`

    Args:
        dcm (sympy.matrices.dense.MutableDenseMatrix):
            Direction cosine matrix transforming vector components from frame A to
            frame B
        diffmap (dict):
            Dictionary of dependent variables and their first derivatives
        diffby (sympy.core.symbol.Symbol):
            Inependent variable to differentiate with respect to. Defaults to t

    Returns:
        sympy.matrices.dense.MutableDenseMatrix:
            Angular velocity in components of frame B

    """

    # s = solve(dcm*difftotalmat(dcm,t,diffmap).T-skew([w1,w2,w3]),(w1,w2,w3))
    # return Matrix([s[w1],s[w2],s[w3]])

    tmp = dcm * difftotalmat(dcm, t, diffmap).T
    return simplify(Matrix([tmp[2, 1], tmp[0, 2], tmp[1, 0]]))


def DCM2axang(DCM):
    r"""Given a direction cosine matrix :math:`{}^\mathcal{B}C^\mathcal{A}` compute
    the axis and angle of the rotation.  Inverse of `calcDCM`.

    Args:
        DCM (sympy.matrices.dense.MutableDenseMatrix):
            Direction cosine matrix transforming vector components from frame A to
            frame B

    Returns:
        tuple:
            n (sympy.matrices.dense.MutableDenseMatrix):
                3x1 matrix representation of the unit vector of the axis of rotation
            th (sympy.core.*):
                Expression for the angle of rotation

    """

    costh = (DCM.trace() - 1) / 2
    sinth = sqrt(1 - costh**2)
    tmp = Matrix([DCM[2, 1] - DCM[1, 2], DCM[0, 2] - DCM[2, 0], DCM[1, 0] - DCM[0, 1]])
    n = tmp / 2 / sinth
    th = atan2(sinth, costh)

    return n, th


def mat2vec(mat, basis="e", hat=True):
    r"""Transform matrix representation of a vector to the vector equation
    for a given basis

    Args:
        mat (sympy.matrices.dense.MutableDenseMatrix):
            3x1 Matrix representation of a vector in components of some frame
        basis (str or iterable):
            If a string, generate unit vector basis for the frame as basis_i (e.g. 'e'
            becomes basis e_1,e_2,e_3). If an iterable of strings (must be of length 3)
            use directly as the basis representation.
        hat (bool):
            Only applies if `basis` input is a string.  If set, basis vectors are
            typeset as bold and hatted (e.g. :math:`\mathbf{\hat{e}}_1`.  If false,
            vectors are only bolded.  Defaults True.
    """

    assert isinstance(basis, str) or (
        hasattr(basis, "__iter__") and len(basis) == 3
    ), "basis input must be a string or iterable of length 3."

    if isinstance(basis, str):
        if hat:
            basis = [r"\hat{\mathbf{" + basis + "}}_" + str(j) for j in range(1, 4)]
        else:
            basis = [r"\mathbf{" + basis + "}_" + str(j) for j in range(1, 4)]

    basissyms = symbols(basis, commutative=False)
    basisvec = Matrix(basissyms)

    return (mat.T * basisvec)[0]


def fancyMat(prefix, shape):
    r"""Create an indexed 2D matrix akin to symarray
    Args:
        prefix (str):
            string representation of symbol to use as the matrix contents
        shape (iterable):
            2-element iterable of the shape of the matrix

    Returns:
        sympy.matrices.dense.MutableDenseMatrix:
            The matrix

    .. note::

        Indexing is 1-based.

    Examples:
        >>> fancyMat('{}^\mathcal{B}C^{\mathcal{A}}',(3,3))
        Matrix([
        [{}^\mathcal{B}C^{\mathcal{A}}_{11}, {}^\mathcal{B}C^{\mathcal{A}}_{12}, {}^\mathcal{B}C^{\mathcal{A}}_{13}],
        [{}^\mathcal{B}C^{\mathcal{A}}_{21}, {}^\mathcal{B}C^{\mathcal{A}}_{22}, {}^\mathcal{B}C^{\mathcal{A}}_{23}],
        [{}^\mathcal{B}C^{\mathcal{A}}_{31}, {}^\mathcal{B}C^{\mathcal{A}}_{32}, {}^\mathcal{B}C^{\mathcal{A}}_{33}]])

    """

    M = []
    for r in range(1, shape[0] + 1):
        row = []
        for c in range(1, shape[1] + 1):
            row.append(prefix + "_{" + str(r) + str(c) + "}")
        M.append(row)
    M = Matrix(symbols(M))

    return M


def fancyVec(prefix, n):
    r"""Create an indexed column matrix akin to symarray

    Args:
        prefix (str):
            string representation of symbol to use as the matrix contents
        n (int):
            Dimension of vector

    Returns:
        sympy.matrices.dense.MutableDenseMatrix:
            The matrix representation of the vector

    .. note::

        Indexing is 1-based.

    Examples:
        >>> fancyVec('a',3)
        Matrix([
        [a_{1}],
        [a_{2}],
        [a_{3}]])

    """

    M = []
    for r in range(1, n + 1):
        M.append(prefix + "_{" + str(r) + "}")
    M = Matrix(symbols(M))

    return M


def EulerAngSet(rots, angs):
    r"""Calculate the equivalent direction cosine matrix for a body Euler Angle set


    Args:
        rots (iterable):
            3-element iterable defining order of rotations of a body Euler angle set.
            For example, a Body-2 3-1-3 rotation would be [3,1,3] and a Body-3 3-2-1
            rotation would be [3,2,1].
        angs (iterable):
            3-elements iterable of symbols or expressions defining the angle of each
            rotation.

    Returns:
        sympy.matrices.dense.MutableDenseMatrix:
            The equivalent direction cosine matrix :math:`{}^\mathcal{B}C^\mathcal{A}`

    """

    assert (
        hasattr(rots, "__iter__") and len(rots) == 3
    ), "rots must be an iterable of length 3."

    assert (
        hasattr(angs, "__iter__") and len(angs) == 3
    ), "v must be an iterable of length 3."

    DCM = eye(3)
    for rot, ang in zip(rots, angs):
        DCM = rotMat(rot, ang) * DCM

    return simplify(DCM)


def EulerLagrange(L, qs, diffmap, diffby=t):
    r"""Apply the Euler-Lagrange equations

    Args:
        L (sympy.core.*):
            Any expression representing the Lagrangian
        qs (iterable):
            The generalized coordinates.  Must be an iterable even if there is only one
            coordinate
        diffmap (dict):
            Dictionary of dependent variables and their first derivatives
        diffby (sympy.core.symbol.Symbol):
            Inependent variable to differentiate with respect to. Defaults to t

    Returns:
        dict:
            Equations of motion

    """

    assert hasattr(qs, "__iter__"), "qs must be an iterable."

    eqs = [simplify(difftotal(L.diff(diffmap[q]), t, diffmap) - L.diff(q)) for q in qs]

    return simplify(solve(eqs, [diffmap[diffmap[q]] for q in qs]))
