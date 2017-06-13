from display import *
from matrix import *
from math import *
from gmath import *


def magnitude(v):
    return math.sqrt(sum([pow(component, 2) for component in v]))

def normalize(v):
    m = magnitude(v)
    return [component / float(m) for component in v]

def subtract_vector(v_one, v_two):
    diff = [0, 0, 0]
    for i in range(3):
        diff[i] = v_one[i] - v_two[i]

    return diff

def scalar_mult_two(v, scalar):
    return [component * scalar for component in v]

def dot_product(v_one, v_two):
    v_one = normalize(v_one)
    v_two = normalize(v_two)
    return v_one[0] * v_two[0] + v_one[1] * v_two[1] + v_one[2] * v_two[2]

def ambient_color(color, ka):
    return color * ka

def diffuse_color(color, kd, normal, light_vector):
    return color * kd * max(0, dot_product(normal, light_vector))

def specular_color(color, ks, normal, light_vector):
    P = scalar_mult_two(normal, dot_product(normal, light_vector))
    R = subtract_vector(scalar_mult_two(P, 2), light_vector)
    view_vector = [0, 0, 1]
    exp = 120

    return color * ks * pow(max(0, dot_product(R, view_vector)), exp)

def illumination(matrix, index, normal, shading, color):
    color = [0, 0, 0]
    
    if len(shading['constants']) > 0:
        for c in shading['constants']:
            ka = [shading['constants'][c]['red'][0], shading['constants'][c]['green'][0], shading['constants'][c]['blue'][0]]
            kd = [shading['constants'][c]['red'][1], shading['constants'][c]['green'][1], shading['constants'][c]['blue'][1]]
            ks = [shading['constants'][c]['red'][2], shading['constants'][c]['green'][2], shading['constants'][c]['blue'][2]]

    for i in range(3):

        location = shading['light']['l1']['location']
        intensity = shading['light']['l1']['color']

        if (len(shading['ambient']) > 0):
            ambient = ambient_color(shading['ambient'][i], ka[i])
        else:
            ambient = ambient_color(intensity[i], ka[i])
            
        color[i] += ambient

        diffuse = diffuse_color(intensity[i], kd[i], normal, location)
        color[i] += diffuse

        specular = specular_color(intensity[i], ks[i], normal, location)
        color[i] += specular

    for i in range(3):
        if color[i] > 255:
            color[i] = 255
        else:
            color[i] = int(round(color[i]))

    return color

def scanline_convert(polygons, i, screen, zbuffer, color, normal, shading):
    if len(shading['light']) > 0 or len(shading['ambient']) > 0:
        color = illumination(polygons, i, normal, shading, color)

    #================================================
    #Done with help from Brian Yang
    x_values = [polygons[i + j][0] for j in range(3)]
    y_values = [polygons[i + j][1] for j in range(3)]
    z_values = [polygons[i + j][2] for j in range(3)]
    pos = [0, 1, 2]

    b_ind = y_values.index(min(y_values))
    t_ind = y_values.index(max(y_values[::-1]))

    pos.remove(b_ind)
    pos.remove(t_ind)
    m_ind = pos[0]
    #================================================

    bot_y = y_values[b_ind]
    top_y = y_values[t_ind]
    mid_y = y_values[m_ind]

    bot_x = x_values[b_ind]
    top_x = x_values[t_ind]
    mid_x = x_values[m_ind]

    bot_z = z_values[b_ind]
    top_z = z_values[t_ind]
    mid_z = z_values[m_ind]

    x0 = x1 = bot_x
    z0 = z1 = bot_z
    delta_x0 = delta_x1 = 0
    delta_z0 = delta_z1 = 0
    
    if top_y != bot_y:
        delta_x0 = float(top_x - bot_x) / (top_y - bot_y)
        delta_z0 = float(top_z - bot_z) / (top_y - bot_y)

    if mid_y != bot_y:
        delta_x1 = float(mid_x - bot_x) / (mid_y - bot_y)
        delta_z1 = float(mid_z - bot_z) / (mid_y - bot_y)
    
    for y in range(int(bot_y), int(mid_y)):
        draw_line(int(x0), int(y), int(z0), int(x1), int(y), int(z1), screen, zbuffer, color)
        x0 += delta_x0
        x1 += delta_x1
        z0 += delta_z0
        z1 += delta_z1

    x1 = mid_x
    z1 = mid_z
    delta_x1 = 0
    delta_z1 = 0
    
    if top_y != mid_y:
        delta_x1 = float(top_x - mid_x) / (top_y - mid_y)
        delta_z1 = float(top_z - mid_z) / (top_y - mid_y)

    for y in range(int(mid_y), int(top_y)):
        draw_line(int(x0), int(y), int(z0), int(x1), int(y), int(z1), screen, zbuffer, color)
        x0 += delta_x0
        x1 += delta_x1
        z0 += delta_z0
        z1 += delta_z1

def add_polygon( polygons, x0, y0, z0, x1, y1, z1, x2, y2, z2 ):
    add_point(polygons, x0, y0, z0);
    add_point(polygons, x1, y1, z1);
    add_point(polygons, x2, y2, z2);

def draw_polygons( matrix, screen, zbuffer, color, shading ):
    if len(matrix) < 2:
        print 'Need at least 3 points to draw'
        return

    point = 0
    while point < len(matrix) - 2:

        normal = calculate_normal(matrix, point)[:]
        #print normal
        if normal[2] > 0:
            polygon_color = [50, 50, 50]
            scanline_convert(matrix, point, screen, zbuffer, polygon_color, normal, shading)

        point+= 3

def add_box( polygons, x, y, z, width, height, depth ):
    x1 = x + width
    y1 = y - height
    z1 = z - depth

    #front
    add_polygon(polygons, x, y, z, x1, y1, z, x1, y, z);
    add_polygon(polygons, x, y, z, x, y1, z, x1, y1, z);

    #back
    add_polygon(polygons, x1, y, z1, x, y1, z1, x, y, z1);
    add_polygon(polygons, x1, y, z1, x1, y1, z1, x, y1, z1);

    #right side
    add_polygon(polygons, x1, y, z, x1, y1, z1, x1, y, z1);
    add_polygon(polygons, x1, y, z, x1, y1, z, x1, y1, z1);
    #left side
    add_polygon(polygons, x, y, z1, x, y1, z, x, y, z);
    add_polygon(polygons, x, y, z1, x, y1, z1, x, y1, z);

    #top
    add_polygon(polygons, x, y, z1, x1, y, z, x1, y, z1);
    add_polygon(polygons, x, y, z1, x, y, z, x1, y, z);
    #bottom
    add_polygon(polygons, x, y1, z, x1, y1, z1, x1, y1, z);
    add_polygon(polygons, x, y1, z, x, y1, z1, x1, y1, z1);

def add_sphere( edges, cx, cy, cz, r, step ):
    points = generate_sphere(cx, cy, cz, r, step)
    num_steps = int(1/step+0.1)

    lat_start = 0
    lat_stop = num_steps
    longt_start = 0
    longt_stop = num_steps

    num_steps+= 1
    for lat in range(lat_start, lat_stop):
        for longt in range(longt_start, longt_stop):

            p0 = lat * (num_steps) + longt
            p1 = p0+1
            p2 = (p1+num_steps) % (num_steps * (num_steps-1))
            p3 = (p0+num_steps) % (num_steps * (num_steps-1))

            if longt != num_steps - 2:
	        add_polygon( edges, points[p0][0],
		             points[p0][1],
		             points[p0][2],
		             points[p1][0],
		             points[p1][1],
		             points[p1][2],
		             points[p2][0],
		             points[p2][1],
		             points[p2][2])
            if longt != 0:
	        add_polygon( edges, points[p0][0],
		             points[p0][1],
		             points[p0][2],
		             points[p2][0],
		             points[p2][1],
		             points[p2][2],
		             points[p3][0],
		             points[p3][1],
		             points[p3][2])

def generate_sphere( cx, cy, cz, r, step ):
    points = []
    num_steps = int(1/step+0.1)

    rot_start = 0
    rot_stop = num_steps
    circ_start = 0
    circ_stop = num_steps

    for rotation in range(rot_start, rot_stop):
        rot = step * rotation
        for circle in range(circ_start, circ_stop+1):
            circ = step * circle

            x = r * math.cos(math.pi * circ) + cx
            y = r * math.sin(math.pi * circ) * math.cos(2*math.pi * rot) + cy
            z = r * math.sin(math.pi * circ) * math.sin(2*math.pi * rot) + cz

            points.append([x, y, z])
            #print 'rotation: %d\tcircle%d'%(rotation, circle)
    return points

def add_torus( edges, cx, cy, cz, r0, r1, step ):
    points = generate_torus(cx, cy, cz, r0, r1, step)
    num_steps = int(1/step+0.1)

    lat_start = 0
    lat_stop = num_steps
    longt_start = 0
    longt_stop = num_steps

    for lat in range(lat_start, lat_stop):
        for longt in range(longt_start, longt_stop):

            p0 = lat * (num_steps) + longt;
            if (longt == num_steps - 1):
	        p1 = p0 - longt;
            else:
	        p1 = p0 + 1;
            p2 = (p1 + num_steps) % (num_steps * num_steps);
            p3 = (p0 + num_steps) % (num_steps * num_steps);

            add_polygon(edges,
                        points[p0][0],
                        points[p0][1],
                        points[p0][2],
                        points[p3][0],
                        points[p3][1],
                        points[p3][2],
                        points[p2][0],
                        points[p2][1],
                        points[p2][2] )
            add_polygon(edges,
                        points[p0][0],
                        points[p0][1],
                        points[p0][2],
                        points[p2][0],
                        points[p2][1],
                        points[p2][2],
                        points[p1][0],
                        points[p1][1],
                        points[p1][2] )

def generate_torus( cx, cy, cz, r0, r1, step ):
    points = []
    num_steps = int(1/step+0.1)

    rot_start = 0
    rot_stop = num_steps
    circ_start = 0
    circ_stop = num_steps

    for rotation in range(rot_start, rot_stop):
        rot = step * rotation
        for circle in range(circ_start, circ_stop):
            circ = step * circle

            x = math.cos(2*math.pi * rot) * (r0 * math.cos(2*math.pi * circ) + r1) + cx;
            y = r0 * math.sin(2*math.pi * circ) + cy;
            z = -1*math.sin(2*math.pi * rot) * (r0 * math.cos(2*math.pi * circ) + r1) + cz;

            points.append([x, y, z])
    return points

def add_circle( points, cx, cy, cz, r, step ):
    x0 = r + cx
    y0 = cy
    t = step

    while t <= 1.00001:
        x1 = r * math.cos(2*math.pi * t) + cx;
        y1 = r * math.sin(2*math.pi * t) + cy;

        add_edge(points, x0, y0, cz, x1, y1, cz)
        x0 = x1
        y0 = y1
        t+= step

def add_curve( points, x0, y0, x1, y1, x2, y2, x3, y3, step, curve_type ):

    xcoefs = generate_curve_coefs(x0, x1, x2, x3, curve_type)[0]
    ycoefs = generate_curve_coefs(y0, y1, y2, y3, curve_type)[0]

    t = step
    while t <= 1.00001:
        x = xcoefs[0] * t*t*t + xcoefs[1] * t*t + xcoefs[2] * t + xcoefs[3]
        y = ycoefs[0] * t*t*t + ycoefs[1] * t*t + ycoefs[2] * t + ycoefs[3]

        add_edge(points, x0, y0, 0, x, y, 0)
        x0 = x
        y0 = y
        t+= step

def draw_lines( matrix, screen, zbuffer, color ):
    if len(matrix) < 2:
        print 'Need at least 2 points to draw'
        return

    point = 0
    while point < len(matrix) - 1:
        draw_line( int(matrix[point][0]),
                   int(matrix[point][1]),
                   matrix[point][2],
                   int(matrix[point+1][0]),
                   int(matrix[point+1][1]),
                   matrix[point+1][2],
                   screen, zbuffer, color)
        point+= 2

def add_edge( matrix, x0, y0, z0, x1, y1, z1 ):
    add_point(matrix, x0, y0, z0)
    add_point(matrix, x1, y1, z1)

def add_point( matrix, x, y, z=0 ):
    matrix.append( [x, y, z, 1] )




def draw_line( x0, y0, z0, x1, y1, z1, screen, zbuffer, color ):

    #swap points if going right -> left
    if x0 > x1:
        xt = x0
        yt = y0
        zt = z0
        x0 = x1
        y0 = y1
        z0 = z1
        x1 = xt
        y1 = yt
        z1 = zt

    x = x0
    y = y0
    z = z0
    A = 2 * (y1 - y0)
    B = -2 * (x1 - x0)
    wide = False
    tall = False

    if ( abs(x1-x0) >= abs(y1 - y0) ): #octants 1/8
        wide = True
        loop_start = x
        loop_end = x1
        dx_east = dx_northeast = 1
        dy_east = 0
        d_east = A
        distance = x1 - x
        dz = (z1 - z0) / abs(float(x1 - x0)) if x1 != x0 else 0
        if ( A > 0 ): #octant 1
            d = A + B/2
            dy_northeast = 1
            d_northeast = A + B
        else: #octant 8
            d = A - B/2
            dy_northeast = -1
            d_northeast = A - B

    else: #octants 2/7
        tall = True
        dx_east = 0
        dx_northeast = 1
        distance = abs(y1 - y)
        dz = (z1 - z0) / abs(float(y1 - y0)) if y1 != y0 else 0
        if ( A > 0 ): #octant 2
            d = A/2 + B
            dy_east = dy_northeast = 1
            d_northeast = A + B
            d_east = B
            loop_start = y
            loop_end = y1
        else: #octant 7
            d = A/2 - B
            dy_east = dy_northeast = -1
            d_northeast = A - B
            d_east = -1 * B
            loop_start = y1
            loop_end = y

    while ( loop_start < loop_end ):
        plot( screen, zbuffer, color, x, y, z )
        if ( (wide and ((A > 0 and d > 0) or (A < 0 and d < 0))) or
             (tall and ((A > 0 and d < 0) or (A < 0 and d > 0 )))):
            x+= dx_northeast
            y+= dy_northeast
            d+= d_northeast
        else:
            x+= dx_east
            y+= dy_east
            d+= d_east

        z += dz
        loop_start+= 1

    plot( screen, zbuffer, color, x, y, z )
