.. _usage:

Usage
============
``sympyhelpers`` is a very thin wrapper around ``sympy`` and so the majority of operations are identical to what they would be if using ``sympy`` directly.  The one exception is differentiation.  ``sympyhelpers`` provides two routines (:py:meth:`~sympyhelpers.sympyhelpers.difftotal` and :py:meth:`~sympyhelpers.sympyhelpers.difftotalmat`) that allow for differentiation using Newton's (dot) notation as opposed to sympy's default Leibniz (ratio) notation.  

Basic Example
----------------
Consider the example of differentiating :math:`\cos(\theta)` with respect to time, where :math:`\theta` is a time-varying quantity.  In regular ``sympy``, this operation would be:

.. code-block:: python

    from sympy import Function,symbols,diff,cos
    t = symbols("t")
    th = Function("theta")(t)
    diff(cos(th),t)

The output would then look like:

.. math::
   \displaystyle - \sin{\left(\theta{\left(t \right)} \right)} \frac{d}{d t} \theta{\left(t \right)}

The equivalent operation with ``sympyhelpers`` is as follows:

.. code-block:: python
   
    from sympyhelpers.sympyhelpers import difftotal, symbols
    th,thd,t = symbols("theta,thetadot,t")
    diffmap = {th:thd}
    difftotal(cos(th),t,diffmap)

.. math::
   \displaystyle - \dot{\theta} \sin{\left(\theta \right)}


:py:meth:`~sympyhelpers.sympyhelpers.difftotalmat` does exactly the same thing, except that it operates on matrices and differentiates the contents of a matrix element by element. 

The utility is that the second version is significantly legible, especially in cases of more complex expressions containing many terms and many first- and second-order derivatives. The drawback, however, is the requirement of defining separate symbol objects for each order of derivative of each variable, and a differentiation map (the ``diffmap`` dictionary in the example, above) that specifies how these are related.  Fortunately, we can automate much of this process. 

Automatically Defining Variables and Differentiation Maps
--------------------------------------------------------------

In the vast majority of cases in classical mechanics, we care about the derivatives of coordinates up through second order.  That is, if we have a generalized coordinate :math:`q`, our expressions will be functions of :math:`q, \dot q, \ddot q` and the relevant differentiation map will be :math:`q \rightarrow \dot q, \dot q \rightarrow \ddot q`. ``sympyhelpers`` provides a utility method ( :py:meth:`~sympyhelpers.sympyhelpers.gendiffvars`) to automatically create both all of the required symbols and the differentiation maps between them. 

:py:meth:`~sympyhelpers.sympyhelpers.gendiffvars` takes one positional input argument, ``syms``, which is a list (or any iterable) of coordinate definitions.  A full coordinate definition is a 3-element tuple (or any other iterable).  The first element of the coordinate definition is the desired variable name.  The second element is the desired coordinate name (that is, the `names` input to :py:func:`~sympy.core.symbol.symbols`).  The third element is the highest order of derivative desired.  If the third input is omitted, then 2 is automatically assumed.  If, instead of providing a tuple, only a string is given for a coordinate definition, then the method assumes that the variable and coordinate name are the same (and that the highest order derivative is 2).  Thus, a coordinate definition of ``'x'`` is identical to ``('x', 'x', 2)``, and a coordinate definition of ``('th', 'theta')`` is identical to ``('th', 'theta', 2)``. 

For each provided coordinate definition, :py:meth:`~sympyhelpers.sympyhelpers.gendiffvars` will generate symbols for all required orders of derivatives along with a differentiation map. The method returns a dictionary of all of the symbol objects and their variable names as well as the complete differentiation map.  For example, consider the case of two coordinates, :math:`\theta, \phi`, where :math:`\dot \theta` is known to be constant.  The setup for this system would be:

.. code-block:: python

   allsyms, diffmap = gendiffvars([('th','theta',1), ('ph', 'phi')])
   locals().update(allsyms)

The update of the `locals` dictionary is optional, but is very convenient, as it allows for direct use of the produced variables in the local namespace. After running this code, the variables ``th, thd, ph, phd, phdd`` will be available for use in the local scope from which it was called, and the ``diffmap`` dictionary will be: ``{theta: thetadot, phi: phidot, phidot: phiddot}``.

 
Vector Representation
--------------------------------

All (Euclidean) vectors are encoded as 3x1 column matrices, which carry along an implicit frame definition.  Given a reference frame :math:`\mathcal I = ({\hat{\mathbf{e}}}_1, {\hat{\mathbf{e}}}_2, {\hat{\mathbf{e}}}_3)` and a vector :math:`\mathbf{r}_{P/O} = x{\hat{\mathbf{e}}}_1 + y{\hat{\mathbf{e}}}_2 + z{\hat{\mathbf{e}}}_3`, we have the equivalency:

  .. math::
    \left[\mathbf{r}_{P/O} \right]_\mathcal{I} = \begin{bmatrix} x \\y\\z\end{bmatrix}_\mathcal{I} \equiv  x{\hat{\mathbf{e}}}_1 + y{\hat{\mathbf{e}}}_2 + z{\hat{\mathbf{e}}}_3

We would define the equivalent ``sympy`` matrix as:

.. code-block:: python

    x, y, z = symbols("x, y, z")
    r_PO_S = Matrix([x, y, z])

In order to convert back from the matrix form to the usual vector component form, ``sympyhelpers`` provides a method :py:meth:`~sympyhelpers.sympyhelpers.mat2vec`.  For the example above, we could run:

.. code-block:: python

    mat2vec(r_PO_S)

To generate :math:`x{\hat{\mathbf{e}}}_1 + y{\hat{\mathbf{e}}}_2 + z{\hat{\mathbf{e}}}_3`.  By default :py:meth:`~sympyhelpers.sympyhelpers.mat2vec` assumes the standard basis :math:`({\hat{\mathbf{e}}}_1, {\hat{\mathbf{e}}}_2, {\hat{\mathbf{e}}}_3)`, but this can be changed via the ``basis`` keyword.  If ``basis`` is set to a string (e.g., ``basis='b'``), then it the method will automatically geenrate an indexed basis set - that is, for ``basis='b'``, the basis will be :math:`({\hat{\mathbf{b}}}_1, {\hat{\mathbf{b}}}_2, {\hat{\mathbf{b}}}_3)`. If you prefer not to put hats on your basis vectors, you can set keyword input ``hat = False``. Alternatively, you can provide a 3-element iterable to the ``basis`` keyword to define your own custom basis vector set. 

``sympyhelpers`` provides two pre-defined basis vector sets:

* ``polarframe``: :math:`\displaystyle  \hat{\mathbf{e}}_r , \hat{\mathbf{e}}_\theta , \hat{\mathbf{e}}_3`
* ``sphericalframe``: :math:`\displaystyle  \hat{\mathbf{e}}_\phi,  \hat{\mathbf{e}}_\theta, \hat{\mathbf{e}}_\rho`

Two additional sets (``polarframe_nohat`` and ``sphericalframe_nohat``) provide the same definitions, without the hats. 
