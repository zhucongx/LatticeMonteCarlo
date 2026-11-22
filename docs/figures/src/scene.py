# Small helpers for the atomistic renderings of the documentation (Blender, Cycles), in the look of the figures of npj Comput. Mater. 11,
# 365 (2025): spheres with a soft highlight, one light from the upper left, strong ambient occlusion, black cell edges, mild perspective.
# Colors as in the line drawings: Al gray, Mg orange, Zn blue; red marks the vacancy and the jumping atom.
# Used by the render_*.py scripts in this folder:   blender -b -P render_xxx.py -- <output folder>
import json
import math

import bpy
from bpy_extras.object_utils import world_to_camera_view
from mathutils import Vector

HEX = {"Al": "AEB3BA", "Mg": "DE8D3A", "Zn": "4A80BD", "red": "C8423B", "blue": "3F78B5", "green": "3E9B5C",
       "yellow": "E0A21B", "dark": "3A3A3A", "gray": "9A9A9A", "white": "FFFFFF", "black": "000000"}
RADIUS = {"Al": 1.00, "Mg": 1.12, "Zn": 0.92, "X": 1.0}      # relative sizes: Mg drawn larger and Zn smaller than Al


def rgb(name):
    """sRGB hex (or a key of HEX) to linear RGB"""
    h = HEX.get(name, name)
    out = []
    for i in (0, 2, 4):
        c = int(h[i:i + 2], 16) / 255
        out.append(c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4)
    return out


def mix(name1, name2, t):
    a, b = HEX.get(name1, name1), HEX.get(name2, name2)
    return "".join(f"{round(int(a[i:i+2], 16) * (1 - t) + int(b[i:i+2], 16) * t):02X}" for i in (0, 2, 4))


def new_scene(width, height, samples=128):
    bpy.ops.wm.read_factory_settings(use_empty=True)
    scene = bpy.context.scene
    scene.render.engine = "CYCLES"
    scene.render.resolution_x, scene.render.resolution_y = width, height
    scene.render.resolution_percentage = 100
    scene.render.film_transparent = True
    scene.render.image_settings.file_format = "PNG"
    scene.render.image_settings.color_mode = "RGBA"
    scene.view_settings.view_transform = "Standard"
    scene.cycles.samples = samples
    scene.cycles.use_denoising = True
    scene.cycles.max_bounces = 8
    scene.cycles.diffuse_bounces = 1                 # little light between the atoms: dark gaps, as with ambient occlusion
    scene.cycles.transparent_max_bounces = 64
    try:
        prefs = bpy.context.preferences.addons["cycles"].preferences
        prefs.compute_device_type = "METAL"
        prefs.get_devices()
        for d in prefs.devices:
            d.use = d.type == "METAL"
        scene.cycles.device = "GPU"
    except Exception as e:  # falls back to the processor
        print("GPU not available:", e)
    world = bpy.data.worlds.new("world")
    scene.world = world
    world.use_nodes = True
    bg = world.node_tree.nodes["Background"]
    bg.inputs["Color"].default_value = (1, 1, 1, 1)
    bg.inputs["Strength"].default_value = 0.42
    return scene


_materials, _meshes = {}, {}


def reset():
    _materials.clear()
    _meshes.clear()


def material(color, alpha=1.0, emission=False):
    key = (color, alpha, emission)
    if key in _materials:
        return _materials[key]
    m = bpy.data.materials.new(f"m{len(_materials)}")
    m.use_nodes = True
    p = m.node_tree.nodes["Principled BSDF"]
    c = rgb(color)
    p.inputs["Base Color"].default_value = (*c, 1)
    p.inputs["Roughness"].default_value = 0.34              # a soft highlight, like the balls of the line drawings
    p.inputs["Specular IOR Level"].default_value = 0.42
    p.inputs["Alpha"].default_value = alpha
    if emission:
        p.inputs["Base Color"].default_value = (0, 0, 0, 1)
        p.inputs["Emission Color"].default_value = (*c, 1)
        p.inputs["Emission Strength"].default_value = 1.0
    _materials[key] = m
    return m


def _outline_material():
    key = "outline"
    if key not in _materials:
        m = bpy.data.materials.new("outline")
        m.use_nodes = True
        nt = m.node_tree
        nt.nodes.clear()
        out = nt.nodes.new("ShaderNodeOutputMaterial")
        mixer = nt.nodes.new("ShaderNodeMixShader")
        tr = nt.nodes.new("ShaderNodeBsdfTransparent")
        em = nt.nodes.new("ShaderNodeEmission")
        em.inputs["Color"].default_value = (0.01, 0.01, 0.01, 1)
        geo = nt.nodes.new("ShaderNodeNewGeometry")
        nt.links.new(geo.outputs["Backfacing"], mixer.inputs[0])
        nt.links.new(tr.outputs[0], mixer.inputs[1])
        nt.links.new(em.outputs[0], mixer.inputs[2])
        nt.links.new(mixer.outputs[0], out.inputs["Surface"])
        _materials[key] = m
    return _materials[key]


def _link(mesh, name):
    ob = bpy.data.objects.new(name, mesh)
    bpy.context.scene.collection.objects.link(ob)
    return ob


def _sphere_mesh(key, mat):
    if key not in _meshes:
        import bmesh
        bm = bmesh.new()
        bmesh.ops.create_uvsphere(bm, u_segments=64, v_segments=32, radius=1.0)
        me = bpy.data.meshes.new("sphere")
        bm.to_mesh(me)
        bm.free()
        for poly in me.polygons:
            poly.use_smooth = True
        me.materials.append(mat)
        _meshes[key] = me
    return _meshes[key]


def sphere(pos, radius, color, alpha=1.0, outline=0.0, emission=False):
    ob = _link(_sphere_mesh(("sphere", color, alpha, emission), material(color, alpha, emission)), "atom")
    ob.location = pos
    ob.scale = (radius,) * 3
    if outline:                                      # a slightly larger sphere of which only the back side is drawn, in black
        h = _link(_sphere_mesh(("outline",), _outline_material()), "outline")
        h.location = pos
        h.scale = (radius + outline,) * 3
        h.visible_diffuse = h.visible_glossy = h.visible_shadow = h.visible_transmission = False
    return ob


def rod(p, q, radius, color, alpha=1.0, emission=False):
    key = ("rod", color, alpha, emission)
    if key not in _meshes:
        import bmesh
        bm = bmesh.new()
        bmesh.ops.create_cone(bm, cap_ends=True, segments=32, radius1=1.0, radius2=1.0, depth=1.0)
        me = bpy.data.meshes.new("rod")
        bm.to_mesh(me)
        bm.free()
        for poly in me.polygons:
            poly.use_smooth = len(poly.vertices) == 4
        me.materials.append(material(color, alpha, emission))
        _meshes[key] = me
    p, q = Vector(p), Vector(q)
    ob = _link(_meshes[key], "rod")
    ob.location = (p + q) / 2
    ob.scale = (radius, radius, (q - p).length)
    ob.rotation_mode = "QUATERNION"
    ob.rotation_quaternion = (q - p).to_track_quat("Z", "Y")
    return ob


def arrow(p, q, radius, color="black", head=0.16, wide=2.8, emission=True):
    """straight arrow from p to q; the head takes the fraction `head` of the length"""
    import bmesh
    p, q = Vector(p), Vector(q)
    m = p.lerp(q, 1 - head)
    rod(p, m, radius, color, 1.0, emission)
    bm = bmesh.new()
    bmesh.ops.create_cone(bm, cap_ends=True, segments=48, radius1=wide * radius, radius2=0.0, depth=(q - m).length)
    me = bpy.data.meshes.new("head")
    bm.to_mesh(me)
    bm.free()
    me.materials.append(material(color, 1.0, emission))
    ob = _link(me, "head")
    ob.location = (m + q) / 2
    ob.rotation_mode = "QUATERNION"
    ob.rotation_quaternion = (q - m).to_track_quat("Z", "Y")


def hull(points, color, alpha):
    """translucent convex hull of the points; returns its edges"""
    import bmesh
    bm = bmesh.new()
    for p in points:
        bm.verts.new(p)
    bmesh.ops.convex_hull(bm, input=bm.verts)
    bmesh.ops.dissolve_limit(bm, angle_limit=0.02, verts=bm.verts, edges=bm.edges)
    edges = [(tuple(e.verts[0].co), tuple(e.verts[1].co)) for e in bm.edges]
    me = bpy.data.meshes.new("hull")
    bm.to_mesh(me)
    bm.free()
    me.materials.append(material(color, alpha))
    ob = _link(me, "hull")
    ob.visible_shadow = False
    return edges


def ring(pos, radius, thickness, color):
    """ring around pos that faces the camera (call it after camera())"""
    import bmesh
    me = bpy.data.meshes.new("ring")
    bm = bmesh.new()
    n, m = 128, 16
    verts = []
    for i in range(n):
        a = 2 * math.pi * i / n
        row = []
        for j in range(m):
            b = 2 * math.pi * j / m
            r = radius + thickness * math.cos(b)
            row.append(bm.verts.new((r * math.cos(a), r * math.sin(a), thickness * math.sin(b))))
        verts.append(row)
    for i in range(n):
        for j in range(m):
            bm.faces.new((verts[i][j], verts[(i + 1) % n][j], verts[(i + 1) % n][(j + 1) % m], verts[i][(j + 1) % m]))
    bm.to_mesh(me)
    bm.free()
    me.materials.append(material(color, 1.0, True))
    ob = _link(me, "ring")
    ob.location = pos
    ob.rotation_mode = "QUATERNION"
    cam = bpy.context.scene.camera
    ob.rotation_quaternion = (cam.location - Vector(pos)).to_track_quat("Z", "Y")
    return ob


def vacancy(pos, radius, color="D0342C"):
    """a vacancy: a white sphere with a red ring that faces the camera"""
    sphere(pos, radius * 0.93, "FFFFFF")
    ring(pos, radius * 0.96, radius * 0.10, color)


def box(lo, hi, radius, color="black", alpha=1.0):
    """edges of the cell, black and unshaded"""
    c = [(x, y, z) for x in (lo[0], hi[0]) for y in (lo[1], hi[1]) for z in (lo[2], hi[2])]
    for i in range(8):
        for j in range(i + 1, 8):
            if sum(abs(a - b) > 1e-9 for a, b in zip(c[i], c[j])) == 1:
                r = rod(c[i], c[j], radius, color, alpha, True)
                r.visible_shadow = r.visible_diffuse = False


def camera(target, azimuth, elevation, scale, perspective=4.5):
    """azimuth 0 looks along +y; scale is the extent that is visible at the target; perspective is the distance in units of scale"""
    az, el = math.radians(azimuth), math.radians(elevation)
    direction = Vector((math.sin(az) * math.cos(el), -math.cos(az) * math.cos(el), math.sin(el)))
    cam = bpy.data.cameras.new("cam")
    distance = perspective * scale
    cam.type = "PERSP"
    cam.sensor_width = 36.0
    cam.angle = 2 * math.atan(0.5 / perspective)
    cam.clip_start, cam.clip_end = 0.1, 10 * distance
    ob = bpy.data.objects.new("cam", cam)
    bpy.context.scene.collection.objects.link(ob)
    ob.location = Vector(target) + direction * distance
    ob.rotation_mode = "QUATERNION"
    ob.rotation_quaternion = direction.to_track_quat("Z", "Y")
    bpy.context.scene.camera = ob
    bpy.context.view_layer.update()
    return ob


def lights(cam, strength=2.9, soft=32.0, offset=(-0.50, 0.62, 1.0)):
    """one sun that shines from the upper left of the camera"""
    light = bpy.data.lights.new("sun", "SUN")
    light.energy = strength
    light.angle = math.radians(soft)
    ob = bpy.data.objects.new("sun", light)
    bpy.context.scene.collection.objects.link(ob)
    ob.rotation_mode = "QUATERNION"
    ob.rotation_quaternion = (cam.rotation_quaternion @ Vector(offset)).to_track_quat("Z", "Y")


def render(path, labels=None):
    """renders to path; labels {name: position} are projected and written next to it as relative coordinates"""
    bpy.context.view_layer.update()
    scene = bpy.context.scene
    if labels:
        out = {}
        for name, pos in labels.items():
            v = world_to_camera_view(scene, scene.camera, Vector(pos))
            out[name] = [round(v.x, 5), round(v.y, 5)]
        with open(path.rsplit(".", 1)[0] + ".json", "w") as f:
            json.dump(out, f)
    scene.render.filepath = path
    bpy.ops.render.render(write_still=True)
    print("rendered", path)


def read_cfg(path):
    """positions (in Angstrom), species and box length of a cubic configuration file"""
    import gzip
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as f:
        lines = f.read().split("\n")
    n = int(lines[0].split("=")[1])
    length = float([l for l in lines if l.startswith("H0(1,1)")][0].split("=")[1].split()[0])
    start = [i for i, l in enumerate(lines) if l.startswith("entry_count")][0] + 1
    body = [l for l in lines[start:] if not l.startswith("auxiliary")]
    species = [body[3 * i + 1].strip() for i in range(n)]
    pos = [tuple(float(v) * length for v in body[3 * i + 2].split()[:3]) for i in range(n)]
    return pos, species, length


def output_folder():
    import sys
    return sys.argv[sys.argv.index("--") + 1:]
