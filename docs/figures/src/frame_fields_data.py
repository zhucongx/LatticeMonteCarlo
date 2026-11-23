# The frame of Ansys that the animation frame_fields shows, for render_frame_turn.py and animate_frame_fields.py: the atoms of its
# clusters with their fields, and the place of each in the drawn box. The frame is step 10 000 000 of the KMC run of docs/README.md,
# analyzed with smallest_cluster_criteria 10 and solvent_bond_criteria 3, file ansys/10000000.xyz.gz.
import gzip
import itertools

CLUSTER = ["8A63A9", "3E9B5C", "D487AC", "E6C343", "2F9AA0", "8B6A4A", "565C66"]     # one color per cluster_id, none is that of a species
BIND, VMAX = ["3F78B5", "F4F4F4", "C8423B"], 0.18      # BIND of style.py, and the range of the color bar in eV
A = 4.046                                              # lattice parameter in Angstrom
NEIGHBORS = [v for v in itertools.product((-1, 0, 1), repeat=3) if sum(c * c for c in v) == 2]


def read_xyz(path):
    """atoms of an extended XYZ file as dictionaries of the fields, and the length of the cubic box"""
    with gzip.open(path, "rt") as f:
        n, comment = int(f.readline()), f.readline()
        fields = comment.split("Properties=")[1].split()[0].split(":")
        length = float(comment.split('Lattice="')[1].split()[0])
        atoms = []
        for _ in range(n):
            words, atom = f.readline().split(), {}
            for name, kind, size in zip(fields[0::3], fields[1::3], map(int, fields[2::3])):
                values = [{"S": str, "I": int, "R": float}[kind](w) for w in words[:size]]
                atom[name], words = (values[0] if size == 1 else values), words[size:]
            atoms.append(atom)
    return atoms, length


def load(path):
    """the atoms that are in a cluster, the length of the box, and the number of contacts in clusters that the faces of the drawn box
    cut. The key xyz of an atom is its place in the drawn box. This box begins, for x, y and z, at the lattice plane at which a face
    cuts the fewest contacts between nearest neighbors of the same cluster and, among these planes, is farthest from all clusters."""
    atoms, L = read_xyz(path)
    planes = round(2 * L / A)                              # lattice planes along an edge, half a lattice parameter apart
    members = [a for a in atoms if a["cluster_id"]]
    for a in members:
        a["site"] = tuple(round(2 * x / A) % planes for x in a["pos"])
    cluster = {a["site"]: a["cluster_id"] for a in members}
    cut = {(k, n): 0 for k in range(3) for n in range(planes)}
    for p, v, k in itertools.product(cluster, NEIGHBORS, range(3)):
        q = tuple((a + b) % planes for a, b in zip(p, v))
        if v[k] == 1 and cluster.get(q) == cluster[p]:
            cut[k, q[k]] += 1                              # the face between the planes of p and q
    faces = []
    for k in range(3):
        room = lambda n: min(min((p[k] - n) % planes, (n - 1 - p[k]) % planes) for p in cluster)
        faces.append(min(range(planes), key=lambda n: (cut[k, n], -room(n))))
    for a in members:
        a["xyz"] = [((n - face) % planes + 0.5) * A / 2 for n, face in zip(a["site"], faces)]
    return members, L, sum(cut[k, n] for k, n in enumerate(faces))


def bound(energy):
    """color of a binding energy, as a hexadecimal number"""
    t = max(-1.0, min(1.0, energy / VMAX))
    a, b = BIND[1], BIND[0 if t < 0 else 2]
    return "".join(f"{round(int(a[i:i + 2], 16) * (1 - abs(t)) + int(b[i:i + 2], 16) * abs(t)):02X}" for i in (0, 2, 4))
