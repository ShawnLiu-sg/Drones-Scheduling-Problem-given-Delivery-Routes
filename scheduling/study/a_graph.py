"""A: unlabelled task chains, graph certificates and branch-and-Benders-cut."""
from collections import deque
from time import perf_counter
import gurobipy as gp
from gurobipy import GRB
from .model import eligible


class Network:
    """Edges are lower differences t[v] >= t[u]+w, optionally activated by a literal."""
    def __init__(self, inst, choices, jobs, follows):
        self.ids = {None: 0}
        self.fixed = []
        self.conditional = {}
        self.ends = []

        def idx(event):
            return self.ids.setdefault(event, len(self.ids))

        def arc(a, b, w, lit=None):
            edge = (idx(a), idx(b), w, lit)
            if lit is None:
                self.fixed.append(edge)
            else:
                self.conditional.setdefault(lit, []).append(edge)

        def eq(a, b, w):
            arc(a, b, w)
            arc(b, a, -w)

        for k, rt in enumerate(inst["routes"]):
            n = len(rt["service"])
            for kind in ("A", "D"):
                eq(None, (k, kind, 0), 0)
            for i in range(1, n):
                eq((k,"D",i-1), (k,"A",i), rt["travel"][i-1])
                arc((k,"A",i), (k,"D",i), rt["service"][i])
            for s, flight in enumerate(rt["sorties"]):
                for kind, pos in (("L",flight["l"]), ("U",flight["b"])):
                    arc((k,"A",pos), (k,kind,s), 0)
                    arc((k,kind,s), (k,"D",pos), 0)
                arc((k,"L",s), (k,"U",s), flight["F"])
                arc((k,"U",s), (k,"L",s), -flight["E"])
                if s:
                    arc((k,"U",s-1), (k,"L",s), 0)
            terminal = (k,"D",n-1)
            self.ends.append(idx(terminal))
            arc(terminal, None, -inst["H"])

        for k, i in jobs:
            tau = inst["routes"][k]["tau"][i]
            arc((k,"R",i), None, 2*tau-inst["H"], ("h",k,i))
            arc((k,"A",i), (k,"R",i), -tau, ("h",k,i))
            arc((k,"R",i), (k,"D",i), tau, ("h",k,i))
        for k, j, i in choices:
            rt = inst["routes"][k]
            p, tau = rt["parcels"][j], rt["tau"][i]
            lit = ("x",k,j,i)
            arc(None, (k,"R",i), p["r"], lit)
            if p["kind"] == "sortie":
                arc((k,"R",i), (k,"L",p["dest"]), tau, lit)
            elif i == p["dest"]:
                arc((k,"R",i), (k,"D",i), tau + rt["service"][i], lit)
        for k,i,l,f in follows:
            arc((k,"R",i), (l,"R",f), 2*inst["routes"][k]["tau"][i],
                ("y",k,i,l,f))
        for v in range(1, len(self.ids)):
            self.fixed.append((0, v, 0, None))

    def evaluate(self, on):
        edges = list(self.fixed)
        for token in sorted(on):
            edges.extend(self.conditional[token])
        n = len(self.ids)
        d, pred = [0]*n, [None]*n
        for _ in range(n):
            last = None
            for e in edges:
                a,b,w,_token = e
                if d[a]+w > d[b]:
                    d[b], pred[b], last = d[a]+w, e, b
            if last is None:
                break
        else:
            # At least one reachable positive cycle. Walk through the tail first.
            v = last
            for _ in range(n):
                v = pred[v][0]
            start, cycle = v, []
            while True:
                e = pred[v]
                cycle.append(e)
                v = e[0]
                if v == start:
                    break
            assert sum(e[2] for e in cycle) > 0
            return None, [{e[3] for e in cycle if e[3] is not None}]
        # Critical path with few activated edges: 0/1 shortest path on tight edges.
        tight = [[] for _ in d]
        for a,b,w,token in edges:
            if d[a]+w == d[b]:
                tight[a].append((b,token))
        cost = [n+1]*n
        cost[0] = 0
        prev, queue = [None]*n, deque([0])
        while queue:
            a = queue.popleft()
            for b, token in tight[a]:
                inc = int(token is not None)
                if cost[a]+inc < cost[b]:
                    cost[b], prev[b] = cost[a]+inc, (a,token)
                    (queue.append if inc else queue.appendleft)(b)
        support = []
        for v in self.ends:
            path = set()
            while v != 0:
                v, token = prev[v]
                if token is not None:
                    path.add(token)
            support.append(path)
        return d, support

    def witness(self, inst, on, d):
        res = []
        for k, rt in enumerate(inst["routes"]):
            assign = [None]*len(rt["parcels"])
            for token in on:
                if token[0] == "x" and token[1] == k:
                    _,_,j,i = token
                    assign[j] = i
            plan = {"assign":assign, "starts":{
                i:d[self.ids[k,"R",i]] for i in sorted(set(assign))}}
            for kind in ("A","D","L","U"):
                n = len(rt["service"] if kind in ("A","D") else rt["sorties"])
                plan[kind] = [d[self.ids[k,kind,i]] for i in range(n)]
            res.append(plan)
        return res


def implied_support(tokens):
    """Both an assignment to e and a selected chain arc incident to e imply h[e]."""
    implied = set()
    for t in tokens:
        if t[0] == "x":
            implied.add(("h",t[1],t[3]))
        elif t[0] == "y":
            implied.update((("h",t[1],t[2]), ("h",t[3],t[4])))
    return tokens-implied


def solve(I, time_limit=30, threads=1):
    """Return a certified bound and best feasible timetable found within the budget."""
    t0 = perf_counter()
    routes = I["routes"]
    if not any(r["parcels"] for r in routes):
        network = Network(I, (), (), ())
        d, _ = network.evaluate(set())
        obj = None if d is None else sum(d[v]-r["base"] for v,r in zip(network.ends,routes))
        return {"status":"infeasible" if d is None else "optimal", "objective":obj,
                "bound":obj, "routes":None if d is None else network.witness(I,set(),d),
                "stats":{"callbacks":0,"cycles":int(d is None),"path_cuts":0,
                         "graph_seconds":perf_counter()-t0,"nodes":0,"master_binary":0}}
    m = gp.Model("A_graph")
    m.Params.OutputFlag = 0
    m.Params.Threads = threads
    m.Params.LazyConstraints = 1
    m.Params.MIPGap = 0
    m.Params.MIPGapAbs = 0
    m.Params.Seed = 0
    choices = [(k,j,i) for k,r in enumerate(routes) for j,p in enumerate(r["parcels"])
               for i in eligible(r,p)]
    jobs = sorted({(k,i) for k,_j,i in choices})
    # Reverse chains on one truck would make the handoff to its earlier stop
    # strictly later than the handoff to its later stop, since tau > 0.
    arcs = [(k,i,l,f) for k,i in jobs for l,f in jobs if k!=l or i<f]
    h = m.addVars(jobs, vtype=GRB.BINARY, name="h")
    x = m.addVars(choices, vtype=GRB.BINARY, name="x")
    y = m.addVars(arcs, vtype=GRB.BINARY, name="y")
    theta = m.addVars(len(routes), lb=0, name="delay")
    # A relaxation in the master: lower bounds on selected departure times
    # along each chain. Actual truck and delivery-sortie timing remains in STN.
    clock = m.addVars(jobs, ub=I["H"], name="chain_start")
    m.setObjective(theta.sum(), GRB.MINIMIZE)
    for k, r in enumerate(routes):
        for j,p in enumerate(r["parcels"]):
            sites = eligible(r,p)
            m.addConstr(x.sum(k,j,"*") == 1)
            # Earliest release + outbound flight + unavoidable truck suffix.
            lbexpr = gp.LinExpr()
            for i in sites:
                suffix = sum(r["travel"][i:]) + sum(r["service"][i+1:])
                own = r["service"][i] if p["kind"]=="truck" and p["dest"]==i else 0
                b = p["r"] + r["tau"][i] + own + suffix - r["base"]
                lbexpr += max(0,b) * x[k,j,i]
            m.addConstr(theta[k] >= lbexpr)
    for k,i in jobs:
        ids = [j for kk,j,ii in choices if kk==k and ii==i]
        rt = routes[k]
        tau = rt["tau"][i]
        earliest_arr = sum(rt["travel"][:i])+sum(rt["service"][:i])
        suffix = sum(rt["travel"][i:])+sum(rt["service"][i+1:])
        m.addConstr(clock[k,i] >= max(0,earliest_arr-tau)*h[k,i])
        m.addConstr(clock[k,i] <= (I["H"]-tau)*h[k,i])
        m.addConstr(theta[k]+rt["base"] >= clock[k,i]+(tau+suffix)*h[k,i])
        m.addConstr(gp.quicksum(routes[k]["parcels"][j]["w"]*x[k,j,i]
                                for j in ids) <= I["Q"]*h[k,i])
        for j in ids:
            m.addConstr(x[k,j,i] <= h[k,i])
            m.addConstr(clock[k,i] >= rt["parcels"][j]["r"]*x[k,j,i])
        m.addConstr(h[k,i] <= x.sum(k,"*",i))
        m.addConstr(y.sum(k,i,"*","*") <= h[k,i])
        m.addConstr(y.sum("*","*",k,i) <= h[k,i])
    m.addConstr(h.sum()-y.sum() <= I["M"])
    for k,i,l,f in arcs:
        tau = routes[k]["tau"][i]
        m.addConstr(clock[l,f] >= clock[k,i]+2*tau-(I["H"]+tau)*(1-y[k,i,l,f]))
    literals = {("x",)+key:var for key,var in x.items()}
    literals.update({("y",)+key:var for key,var in y.items()})
    literals.update({("h",)+key:var for key,var in h.items()})
    keys = list(literals)
    vlist = list(literals.values())
    network = Network(I, choices, jobs, arcs)
    warm = I.get("warm")
    best = [None,None]
    if warm is not None:
        best[:] = [sum(p["D"][-1]-r["base"] for p,r in zip(warm,routes)),warm]
        m.addConstr(theta.sum() <= best[0])
        task_list = sorted((start,k,int(i)) for k,p in enumerate(warm)
                           for i,start in p["starts"].items())
        tails, used_arcs = [], set()
        for t,k,i in task_list:
            a = next((a for a,(free,_) in enumerate(tails) if free<=t),len(tails))
            end = t+2*routes[k]["tau"][i]
            if a==len(tails):
                tails.append((end,(k,i)))
            else:
                used_arcs.add(tails[a][1]+(k,i))
                tails[a] = (end,(k,i))
        for (k,j,i),v in x.items():
            v.Start = int(warm[k]["assign"][j]==i)
        for (k,i),v in h.items():
            v.Start = int(i in warm[k]["assign"])
        for key,v in y.items():
            v.Start = int(key in used_arcs)
        for k,v in theta.items():
            v.Start = warm[k]["D"][-1]-routes[k]["base"]
    counts = {"callbacks":0,"cycles":0,"path_cuts":0,"graph_seconds":0.0}
    errors = []

    def callback(model, where):
        if where != GRB.Callback.MIPSOL:
            return
        try:
            counts["callbacks"] += 1
            on = {k for k,val in zip(keys,model.cbGetSolution(vlist)) if val>.5}
            g0 = perf_counter()
            d,support = network.evaluate(on)
            counts["graph_seconds"] += perf_counter()-g0
            if d is None:
                S = implied_support(support[0])
                model.cbLazy(gp.quicksum(1-literals[i] for i in sorted(S)) >= 1)
                counts["cycles"] += 1
                return
            delays = [d[v]-r["base"] for v,r in zip(network.ends,routes)]
            obj = sum(delays)
            if best[0] is None or obj < best[0]:
                best[:] = [obj, network.witness(I,on,d)]
            for k,v in enumerate(model.cbGetSolution(list(theta.values()))):
                if v+1e-6 < delays[k]:
                    S = implied_support(support[k])
                    diff = gp.quicksum(1-literals[i] for i in sorted(S))
                    model.cbLazy(theta[k] >= delays[k]*(1-diff))
                    counts["path_cuts"] += 1
        except Exception as error:
            errors.append(error)
            model.terminate()

    m.Params.TimeLimit = max(0.0, time_limit-(perf_counter()-t0))
    m.optimize(callback)
    if errors:
        m.dispose()
        raise errors[0]
    status = ("infeasible" if m.Status==GRB.INFEASIBLE else
              "optimal" if m.Status==GRB.OPTIMAL else "limit")
    bound = getattr(m,"ObjBound",-GRB.INFINITY)
    if status=="infeasible" and best[0] is not None:
        raise RuntimeError("Graph master infeasible despite a feasible witness")
    if status=="optimal" and (best[0] is None or abs(best[0]-m.ObjVal)>1e-5):
        raise RuntimeError("Graph solution and final master objective disagree")
    counts.update(nodes=m.NodeCount, master_binary=len(keys), solver_status=m.Status)
    m.dispose()
    return {"status": status, "objective":best[0],
            "bound": None if status=="infeasible" else max(0,bound) if abs(bound)<1e90 else 0,
            "routes":best[1], "stats":counts}
