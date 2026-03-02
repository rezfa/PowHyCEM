# ─── Algorithm & Convergence Parameters ───────────────────────────────────────

const INTEGER_MODE  = Ref(false)

LB        = -Inf
UB        =  Inf
k         = 1
max_iter  = 5000
tolerence = 1e-2
gap       = (UB - LB) / max(1e-6, abs(LB))

const ALPHA        = 0.5
const LEVELSET_CON = Ref{Union{ConstraintRef, Nothing}}(nothing)
φ = 0
