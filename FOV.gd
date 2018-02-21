extends Node2D
	
func _ready():
	run_tests()
	pass

"""
@author: Alex Hall

Computes a 2D Visibility Polygon using Shadowcasting. 

Algorithm from: https://www.redblobgames.com/articles/visibility/
Additional Implementation Ideas from: https://github.com/trylock/visibility

Notes: 
* Not optimized for performance. It is reasonable for low geometry counts, but may not scale well.
* Favors readability over concision
* Cannot handle intersecting polygons. 
	* This could be added by computing segment intersections and dividing into smaller segments at intersection points
* FOV polygon must be closed
	* This can be ensured by adding a bounding box around the scene. 
"""
func get_fov_from_polygons(polygons, origin):
	var segments = polygon_collection_to_segments(polygons, origin)
	var fov = get_fov_from_segments(segments, origin)
	return fov
				
func polygon_collection_to_segments(polygons, origin):
	var fov_segments = []
	for p in polygons:
		for segment in polygon_to_segments(p, origin):
			fov_segments.append(segment)
	return fov_segments
	
func polygon_to_segments(poly, origin):
	#need to pass in origin to filter segments colinear with origin. 
	var segments = []
	for i in range(len(poly)):
		var new_seg = LineSegment.new(poly[i], poly[i-1])
		if valid_segment(new_seg, origin):
			segments.append(new_seg)
	return segments
	
func valid_segment(segment, origin):
	var oab = Comparisons.new(origin).compute_orientation(segment.a, segment.b, origin)
	if oab == ORIENTATION.colinear:
		pass
	return oab != ORIENTATION.colinear
	
func get_fov_from_segments(segments, origin):
	var data = store_segments_as_events(segments, origin, [], OrderedSet.new(origin))
	data['events'] = sort_events_by_angle(origin, data['events'])
	var vertices = find_visibility_polygon(data['events'], data['state'], origin)
	var polygon = remove_colinear_points(vertices)
	return polygon
	
func store_segments_as_events(segment_collection, origin, events=[], state=OrderedSet.new(Vector2(0,0))):
	for segment in segment_collection:
		#sort segments and add them as events
		# skip line segments colinear with point
		var reverse_segment = LineSegment.new(segment.b, segment.a)
		var a = segment.a
		var b = segment.b
		var pab = Comparisons.compute_orientation(origin, a, b)
		if pab == ORIENTATION.colinear:
			print('skipping')
			continue
		elif pab == ORIENTATION.right_turn:
			events.append(VisibilityEvent.new(EVENT_TYPE.start_vertex, segment))
			events.append(VisibilityEvent.new(EVENT_TYPE.end_vertex, reverse_segment))
		else:
			events.append(VisibilityEvent.new(EVENT_TYPE.start_vertex, reverse_segment))
			events.append(VisibilityEvent.new(EVENT_TYPE.end_vertex, segment))
		
		# initialize state by adding line segments that are intersected
		# by vertical ray
		if a.x < b.x:
			var tmp = a
			a = b
			b = tmp

		var abp = Comparisons.compute_orientation(a,b, origin)
		var b_on_south_ray = Comparisons.approx_equal_float(b.x, origin.x)
		
		var south_ray_intersects_segment = a.x > origin.x and origin.x > b.x
		
		var is_clockwise_edge = abp == ORIENTATION.right_turn
		if is_clockwise_edge and (b_on_south_ray or south_ray_intersects_segment):
			state.insert(segment)
		
	var data = {
		'events': events,
		'state': state,
	}	
	return data 
	
func sort_events_by_angle(origin, events):
	events.sort_custom(Comparisons.new(origin), "compare_event")
	return events	
		
func find_visibility_polygon(events, state, origin):
	# ray shoots down 
	# recall: positive Y points DOWN
	var vertices = []
	for event in events:
		if event.type == EVENT_TYPE.end_vertex:
			var s = event.segment
			var rs = reverse_segment(s)
			state.erase(event.segment)
			state.erase(reverse_segment(event.segment))
		if state.empty():
			vertices.append(event.point())
		elif Comparisons.new(origin).compare_segment(event.segment, state.first()):
			# nearest line segment has changed
#			# compute the intersection point with this segment
			var ray = Ray.new(origin, event.point() - origin)
			var nearest_segment = state.first()
			var intersection_data = ray_intersects(ray, nearest_segment)
			var intersection = intersection_data['at']
			var intersects = intersection_data['exists']
			
			#ray intersects line segment L iff L is in state
			assert intersects
			if event.type == EVENT_TYPE.start_vertex:
				vertices.append(intersection)
				vertices.append(event.point())
			else:
				vertices.append(event.point())
				vertices.append(intersection)
		
		if event.type == EVENT_TYPE.start_vertex:
			state.insert(event.segment)
	return vertices
	
func remove_colinear_points(vertices):
	var top = 0
	for idx in range(len(vertices)):
		var prev = top - 1
		var next = (idx + 1) % len(vertices)
		var pvtx = vertices[prev]
		var cvtx = vertices[idx]
		var nvtx = vertices[next]
		var orientation = Comparisons.compute_orientation(pvtx, cvtx, nvtx)
		if orientation != ORIENTATION.colinear:
			vertices[top] = vertices[idx]
			top += 1
	
	var new_verts = []
	for i in range(top):
		new_verts.append(vertices[i])
	
	return new_verts


				
				
# *********************************************************************************************************
# Support Functions
# *********************************************************************************************************
	
func normal(vec2):
	#this is to parallel code from https://github.com/trylock/visibility/blob/master/tests/vector2_test.cpp
	#refactor once working.
	return (-vec2).tangent()

func reverse_segment(line_segment):
	return LineSegment.new(line_segment.b, line_segment.a)

static func ray_intersects(ray, segment):
	var ao = ray.origin - segment.a
	var ab = segment.b - segment.a
	var det = Comparisons.cross(ab, ray.direction)
	var false_result = {
		'exists':false,
		'at':null
	}
	var output = {
		'exists':true,
		'at':null,
	}
	if Comparisons.approx_equal_float(det,0.0):
		var abo = Comparisons.compute_orientation(segment.a, segment.b, ray.origin)
		if abo != ORIENTATION.colinear:
			return false_result
		var dist_a = ao.dot(ray.direction)
		var dist_b = (ray.origin-segment.b).dot(ray.direction)
		if dist_a > 0.0 and dist_b > 0.0:
			return false_result
		elif ((dist_a > 0.0) != (dist_b >0.0)):
			output['at'] = ray.origin
		elif dist_a > dist_b:
			output['at'] = segment.a
		else:
			output['at'] = segment.b
	else:	
		var u = Comparisons.cross(ao, ray.direction) / det
		if u < 0.0 or 1.0 < u:
			return false_result
		
		var t = -1 * Comparisons.cross(ab, ao) / det
		output['at'] = ray.origin + t * ray.direction
		output['exists'] = Comparisons.approx_equal_float(t,0.0) or t > 0.0
	return output
	
	
# *********************************************************************************************************
# Support Classes
# *********************************************************************************************************
const FLOAT_EPSILON = 0.00001

enum ORIENTATION {
        right_turn = 1,
        left_turn = -1,
        colinear = 0
    };
	
enum EVENT_TYPE{start_vertex=0, 
				end_vertex=1}	
	
class VisibilityEvent:

	var type
	var segment
	
	func _init(type=null, segment=null):
		self.type = type
		self.segment = segment
		
	func point():
		return self.segment.a	
		
class LineSegment:
	export var a = Vector2()
	export var b = Vector2()
	func _init(a, b):
		self.a = a
		self.b = b
		
class Ray:
	export var origin = Vector2()
	export var direction = Vector2()
	func _init(origin, direction):
		self.origin = origin
		self.direction = direction		


class Comparisons:
	var origin
	func _init(origin = Vector2(16,16)):
		self.origin = origin
		
	static func approx_equal_float(a, b, epsilon = FLOAT_EPSILON):
	    return abs(a - b) <= epsilon

	static func cross(a, b):
		return a.x * b.y - a.y * b.x
	
	static func approx_equal_vector(a,b,epsilon=FLOAT_EPSILON):
		var same_x = approx_equal_float(a.x, b.x, epsilon)
		var same_y = approx_equal_float(a.y, b.y, epsilon)
		return same_x and same_y 
		
	func compare_segment(a,b):
		return line_segment_is_closer(a,b, self.origin) 
		
	static func compute_orientation(a, b, c):
		var det = cross(b - a, c - a)
		var pos = 1 if 0.0 < det else 0.0
		var neg = 1 if det < 0.0 else 0.0
		return pos - neg

		
	static func line_segment_is_closer(x,y, origin):
		var a = x.a
		var b = x.b
		var c = y.a
		var d = y.b
		assert not compute_orientation(origin, a, b) == ORIENTATION.colinear
		assert not compute_orientation(origin, c, d) == ORIENTATION.colinear
		
		#sort the endpoints so that if there are common endpoints, they will be a and c
		if approx_equal_vector(b,c) or approx_equal_vector(b,d):
			var tmp = a
			a = b
			b = tmp
		if approx_equal_vector(a,d):
			var tmp = c
			c = d
			d = tmp
		
		#cases with common endpoints
		if approx_equal_vector(a, c):
			var oad = compute_orientation(origin, a, d)
			var oab = compute_orientation(origin, a, b)
			var a_is_between_d_and_b = oad != oab
			if approx_equal_vector(b, d) or a_is_between_d_and_b:
				return false
				
			var ab_is_between_origin_and_d = compute_orientation(a, b, d) != compute_orientation(a,b,origin)
			return ab_is_between_origin_and_d
		
		#cases without common endpoints
		var cda = compute_orientation(c,d,a)
		var cdb = compute_orientation(c,d,b)
		var segments_are_colinear_and_parallel_to_origin = cdb == ORIENTATION.colinear and cda == ORIENTATION.colinear
		if segments_are_colinear_and_parallel_to_origin:
			var a_is_closer_to_origin = origin.distance_to_squared(a) < origin.distance_to_squared(c)
			return 	a_is_closer_to_origin	
			
		var cd_and_ab_will_never_intersect = cda == cdb
		var line_from_cd_touches_enpdoint_a = cda == ORIENTATION.colinear
		var line_from_cd_touches_endpoint_b = cdb == ORIENTATION.colinear
		var cd_only_touches_an_endpoint = line_from_cd_touches_endpoint_b or line_from_cd_touches_enpdoint_a
		var no_intersection = cd_and_ab_will_never_intersect or cd_only_touches_an_endpoint
		if no_intersection:
			var cdo = compute_orientation(c, d, origin)
			return cdo == cda or cdo == cdb
		else: #extending cd would intersect ab
			var abo = compute_orientation(a, b, origin)
			var abc = compute_orientation(a,b,c)
			var ab_intersects_cd = abo != abc
			return ab_intersects_cd
			
				
	func compare_event(a,b):
		if approx_equal_vector(a.point(), b.point()):
			return (a.type == EVENT_TYPE.end_vertex and b.type == EVENT_TYPE.start_vertex)
		return is_clockwise(a.point(), b.point(), origin)
		
	static func is_counter_clockwise(a, b, origin):
		# compare angles clockwise starting at the positive y axis
		var a_is_left = a.x < origin.x
		var b_is_left = b.x < origin.x
		if a_is_left != b_is_left:
			return b_is_left
			
		if approx_equal_float(a.x, origin.x) and approx_equal_float(b.x, origin.x):
			if not a.y < origin.y or not b.y < origin.y:
				return b.y < a.y
			return a.y < b.y
		
		var oa = a - origin
		var ob = b - origin
		var det = cross(oa, ob)
		if approx_equal_float(det, 0.0):
			return oa.length_squared() < ob.length_squared()
		return det < 0
		
	static func is_clockwise(a,b,origin):
		var a_is_right = a.x > origin.x
		var b_is_right = b.x > origin.x
		if a_is_right != b_is_right:
			return b_is_right
		if approx_equal_float(a.x, b.x) and approx_equal_float(b.x, origin.x):
			if not a.y > origin.y or not b.y > origin.y:
				return b.y < a.y
			return a.y < b.y
		var oa = a - origin
		var ob = b - origin
		var det = cross(oa, ob)
		if approx_equal_float(det, 0.0):
			return oa.length_squared() > ob.length_squared()
		return det > 0.0 
		
			
class OrderedSet:
	var data = []
	var origin = null
	
	func _init(origin):
		self.origin = origin
		
	func empty():
		return data.empty()
		
	func erase(object):
		var match_idx = find(object)
		if match_idx != null:
				self.data.remove(match_idx)
					
	func has(object):
		return find(object) != false

	func object_matches(x,y):
		return x.a == y.a and x.b == y.b 
	
	func length():
		return len(data)
		
	func find(object):
		var best_idx = bsearch(object, self.data)
		if len(self.data) >  best_idx:
			var obj = self.data[best_idx]
			if 	object_matches(object, obj):
				return best_idx 
		return null
	
	func first():
		return data[0] if not self.empty() else null
		
	func middle(list):
		var mid_idx = floor(len(list) / 2)
		return mid_idx
		
	func bsearch(segment, dataset):
		return dataset.bsearch_custom(segment, Comparisons.new(origin), "compare_segment")
			

	func insert(segment):
		var best_idx = bsearch(segment, self.data)
		if  len(self.data) <= best_idx or not object_matches(segment, self.data[best_idx]):
			self.data.insert(best_idx, segment)


			

			
	
# *********************************************************************************************************
# Tests
# *********************************************************************************************************

func run_tests():
	vector_tests()
	primitive_tests()
	visibility_tests()
	print('FOV tests complete!')


func vector_tests():
	test_create_vector()
	test_add_vectors()
	test_subtract_vectors()
	test_multiply_vectors()
	test_divide_vectors()
	test_negate_vector()
	test_compare_two_vectors()
	test_dot_product()
	test_squared_length()
	test_square_of_euclidean_distance()
	test_2d_normal_vector()
	test_determinant()
	test_normalize_float_vec()
	test_is_counter_clockwise()

func primitive_tests():
	test_ordered_set()
	test_calculate_orientation_of_3_points_in_plane()
	test_calculate_intersection()
	test_line_segment_is_closer_works()
	
func visibility_tests():
	test_compare_line_segments_with_no_common_endpoints()
	test_compare_line_segments_with_common_endpoints()
	test_compare_angle_with_two_points_in_general_position()
	test_compare_angle_with_two_points_if_colinear_with_origin()
	test_sort_events_by_angle()
	test_store_segments_as_events()
	test_storing_segments_as_events()
	test_calculate_visibility_polygon_with_no_line_segments()
	test_calculate_visibility_polygon_with_no_obstacle_only_boundary()
	test_calculate_visibility_polygon_with_a_polyline_obstacle()
	test_calculate_visibility_polygon_with_a_convex_polygon_obstacle()
	test_calculate_visibility_polygon_with_a_concave_polygon_obstacle()

func test_ordered_set():
	var a = LineSegment.new(Vector2(1,1), Vector2(1,0))
	var b = LineSegment.new(Vector2(1,2), Vector2(2,0))
	var c = LineSegment.new(Vector2(1,3), Vector2(3,0))
	var d = LineSegment.new(Vector2(0,4), Vector2(-1,4))
	var set = OrderedSet.new(Vector2(0,0))
	assert set.length() == 0
	set.insert(a)
	assert set.length() == 1
	assert set.find(a) == 0
	set.erase(b)
	assert set.length() == 1
	assert set.find(a) == 0
	set.insert(a)
	assert set.length() == 1
	set.insert(c) 
	assert set.find(c) == 1
	assert set.find(a) == 0
	assert set.length() == 2
	
	set.insert(b)
	assert set.length() == 3
	assert set.find(a) == 0
	assert set.find(b) == 1
	assert set.find(c) == 2
	
	set.erase(b)
	assert set.length() == 2
	assert set.find(a) == 0
	assert set.find(c) == 1
	
	
	
func test_create_vector():
	var result = Vector2(1, 1)
	assert result.x == 1
	assert result.y == 1
	
	var other = Vector2(3, 4)
	assert other.x == 3
	assert other.y == 4
	
func test_add_vectors():
	var a = Vector2(1,2)
	var b = Vector2(3,4)
	
	assert (a + b) == Vector2(4,6)
	assert (a + b).y == 6
	
	a += b
	assert a.x == 4
	assert a.y == 6
	
func test_subtract_vectors():
	var a = Vector2(1,2)
	var b = Vector2(3,4)
	
	assert (a - b) == Vector2(-2,-2)
	a -= b
	assert a.x == -2
	assert a.y == -2
	
func test_multiply_vectors():
	var a = Vector2(1,2)
	assert (a * 3).x == 3
	assert (a * 3).y == 6
	
	a *= 3
	assert a.x == 3
	assert a.y == 6
	
func test_divide_vectors():
	var a = Vector2(2, 8)
	assert (a / 2).x == 1
	assert (a / 2).y == 4
	
	a /= 2 
	assert a.x == 1
	assert a.y == 4
	
func test_negate_vector():
	var a = Vector2(2, 8)
	assert (-a).x == -2
	assert (-a).y == -8
	
func test_compare_two_vectors():
	var a = Vector2(1, 2)
	var b = Vector2(1, 2)
	assert a == b
	
	assert not (a != b)
	
	b.x = 0
	assert not (a == b)
	assert a != b
	
func test_dot_product():
	assert Vector2(1,2).dot(Vector2(3,4)) == 11
	assert Vector2(1,2).dot(Vector2(0,0)) == 0

func test_squared_length():
	assert Vector2(3,4).length_squared() == 25
	assert Vector2(0,0).length_squared() == 0
	
func test_square_of_euclidean_distance():
	assert Vector2(3,4).distance_squared_to(Vector2(0,1)) == 18
	assert Vector2(3,4).distance_squared_to(Vector2(3,4)) == 0 
	
		
func test_2d_normal_vector():
	var a = Vector2(3,4)
	var ortho = normal(a)
	assert ortho.x == -4
	assert ortho.y == 3


	
func test_determinant():
	var a = Vector2(3,4)
	var b = Vector2(1,2)
	var det = Comparisons.cross(a, b)
	assert det == 2

func test_normalize_float_vec():
	var a = Vector2(3,4)
	var normalized = a.normalized()
	assert normalized.length_squared() == 1
	
	var zero = Vector2(0,0)
	normalized = zero.normalized()
	assert normalized.x == 0
	assert normalized.y == 0
	

	
func test_calculate_orientation_of_3_points_in_plane():
	assert Comparisons.compute_orientation(Vector2(0,0), Vector2(1,0), Vector2(2,1)) == ORIENTATION.right_turn
	assert Comparisons.compute_orientation(Vector2(0,0), Vector2(1,0), Vector2(2,-1)) == ORIENTATION.left_turn
	assert Comparisons.compute_orientation(Vector2(0,0), Vector2(1,0), Vector2(2,0)) == ORIENTATION.colinear
	assert Comparisons.compute_orientation(Vector2(0,0), Vector2(0,0), Vector2(4,5)) == ORIENTATION.colinear
	assert Comparisons.compute_orientation(Vector2(0,0), Vector2(0,0), Vector2(0,0)) == ORIENTATION.colinear
			
func test_calculate_intersection():
	var point = Vector2()
	var test_ray = Ray.new(Vector2(0,0), Vector2(1,0))
	var intersection = ray_intersects(test_ray,LineSegment.new(Vector2(-1,1), Vector2(-1,-1)))
	assert not intersection['exists']
	
	intersection = ray_intersects(test_ray, LineSegment.new(Vector2(-.001,1), Vector2(-.001,-1)))
	assert not intersection['exists']
	
	intersection = ray_intersects(test_ray, LineSegment.new(Vector2(-2,0), Vector2(-1,0)))
	assert not intersection['exists']
	intersection = ray_intersects(test_ray,LineSegment.new(Vector2(0,1), Vector2(0,-1)))
	assert intersection['exists']
	point = intersection['at']
	assert point.x == 0
	assert point.y == 0
	
	intersection = ray_intersects(test_ray, LineSegment.new(Vector2(-1, 0), Vector2(0,0)))
	assert intersection['exists']
	point = intersection['at']
	assert point.x == 0
	assert point.y == 0
	
	intersection = ray_intersects(test_ray, LineSegment.new(Vector2(0, 0), Vector2(-1,0)))
	assert intersection['exists']
	point = intersection['at']
	assert point.x == 0
	assert point.y == 0
	
	intersection = ray_intersects(test_ray, LineSegment.new(Vector2(2, 1), Vector2(2,-1)))
	assert intersection['exists']
	point = intersection['at']
	assert point.x == 2
	assert point.y == 0
	
	intersection = ray_intersects(test_ray, LineSegment.new(Vector2(2,0), Vector2(3,0)))
	assert intersection['exists']
	point = intersection['at']
	assert point.x == 2
	assert point.y == 0
	
	intersection = ray_intersects(test_ray, LineSegment.new(Vector2(3,0), Vector2(2,0)))
	assert intersection['exists']
	point = intersection['at']
	assert point.x == 2
	assert point.y == 0
	
	intersection = ray_intersects(test_ray, LineSegment.new(Vector2(1,0), Vector2(1,-1)))
	assert intersection['exists']
	point = intersection['at']
	assert point.x == 1
	assert point.y == 0
	
	intersection = ray_intersects(test_ray, LineSegment.new(Vector2(1,0), Vector2(0,-1)))
	assert intersection['exists']
	point = intersection['at']
	assert point.x == 1
	assert point.y == 0
	
	intersection = ray_intersects(test_ray, LineSegment.new(Vector2(1,0), Vector2(0,1)))
	assert intersection['exists']
	point = intersection['at']
	assert point.x == 1
	assert point.y == 0

	intersection = ray_intersects(
		Ray.new(Vector2(0,0), Vector2(250,-250)),
		LineSegment.new(Vector2(250,-250), Vector2(-250,-250))
	)
	assert intersection['exists']
	point = intersection['at']
	assert point.x == 250
	assert point.y == -250

	

		
func test_line_segment_is_closer_works():
	var origin = Vector2(0,0)
	var a = Vector2(0,0.1)
	var b = Vector2(1,0)
	var c = Vector2(-1,0)
	
	var d = Vector2(0.5,0)
	var e = Vector2(0, 0.5)
	var f = Vector2(0, 1)
	
	var ab = LineSegment.new(a,b)
	var ba = LineSegment.new(b,a)
	var ac = LineSegment.new(a,c)
	
	var de = LineSegment.new(d,e)
	var bf = LineSegment.new(b,f)
	
	assert not Comparisons.line_segment_is_closer(ab, ba,origin)
	assert not Comparisons.line_segment_is_closer(ab, ac,origin)
	assert not Comparisons.line_segment_is_closer(ba, ac,origin)
	assert Comparisons.line_segment_is_closer(de, bf,origin)
	assert not Comparisons.line_segment_is_closer(bf, de,origin)
			
	
	
func test_line_segment_is_closer(a, b, c, d):
	var ab = LineSegment.new(a,b)
	var ba = LineSegment.new(b,a)
	var cd = LineSegment.new(c,d)
	var dc = LineSegment.new(d,c)
	var origin = Vector2(0,0)
	
	assert Comparisons.line_segment_is_closer(ab, cd, origin)
	assert Comparisons.line_segment_is_closer(ba, cd, origin)
	assert Comparisons.line_segment_is_closer(ab, dc, origin)
	assert Comparisons.line_segment_is_closer(ba, dc, origin)
	
	assert not Comparisons.line_segment_is_closer(cd, ab, origin)
	assert not Comparisons.line_segment_is_closer(dc, ab, origin)
	assert not Comparisons.line_segment_is_closer(cd, ba, origin)
	assert not Comparisons.line_segment_is_closer(dc, ba, origin)
	
func test_line_segments_are_equal(a,b,c,d):
	var ab = LineSegment.new(a,b)
	var ba = LineSegment.new(a,b)
	var cd = LineSegment.new(c,d)
	var dc = LineSegment.new(d,c)
	var origin = Vector2(0,0)
	assert not Comparisons.line_segment_is_closer(ab, cd, origin)
	assert not Comparisons.line_segment_is_closer(ba, cd, origin)
	assert not Comparisons.line_segment_is_closer(ab, dc, origin)
	assert not Comparisons.line_segment_is_closer(ba, dc, origin)
	
	assert not Comparisons.line_segment_is_closer(cd, ab, origin)
	assert not Comparisons.line_segment_is_closer(dc, ab, origin)
	assert not Comparisons.line_segment_is_closer(cd, ba, origin)
	assert not Comparisons.line_segment_is_closer(dc, ba, origin)

func test_compare_line_segments_with_no_common_endpoints():
	test_line_segment_is_closer(Vector2(1,1), Vector2(1,-1), Vector2(2,1), Vector2(2,-1))
	test_line_segment_is_closer(Vector2(1,1), Vector2(1,-1), Vector2(2,2), Vector2(2,3))
	
func test_compare_line_segments_with_common_endpoints():
	test_line_segments_are_equal(Vector2(1,1), Vector2(1,0), Vector2(1,0), Vector2(1,-1))
	test_line_segments_are_equal(Vector2(1,1), Vector2(1,0), Vector2(1,0), Vector2(1,1))
	
	test_line_segment_is_closer(Vector2(2,0), Vector2(1,1), Vector2(2,1), Vector2(2,0))
	test_line_segment_is_closer(Vector2(2,1), Vector2(2,0), Vector2(2,0), Vector2(3,2))
	
	
func test_compare_angle_with_two_points_in_general_position():
	var origin = Vector2(0,0)
	assert Comparisons.is_counter_clockwise(Vector2(0,1), Vector2(1,1), origin)
	assert not Comparisons.is_counter_clockwise(Vector2(1,1), Vector2(0,1), origin)
	
	assert Comparisons.is_counter_clockwise(Vector2(1,1), Vector2(1,-1), origin)
	assert not Comparisons.is_counter_clockwise(Vector2(1,-1),Vector2(1,1), origin)
	
	assert Comparisons.is_counter_clockwise(Vector2(1,0), Vector2(-1,-1), origin)
	assert not Comparisons.is_counter_clockwise(Vector2(-1,-1), Vector2(1,0), origin)
	
	assert Comparisons.is_counter_clockwise(Vector2(0,1), Vector2(0,-1), origin)
	assert not Comparisons.is_counter_clockwise(Vector2(0,-1), Vector2(0,1), origin)
	
func test_compare_angle_with_two_points_if_colinear_with_origin():
	var origin = Vector2(0,0)
	assert Comparisons.is_counter_clockwise(Vector2(1,0), Vector2(2,0), origin)
	assert not Comparisons.is_counter_clockwise(Vector2(2,0), Vector2(1,0), origin)
	
	assert not Comparisons.is_counter_clockwise(Vector2(1,0), Vector2(1,0), origin)
	assert not Comparisons.is_counter_clockwise(Vector2(0,0), Vector2(0,0), origin)
	


func test_is_counter_clockwise():
	var origin = Vector2(0,0)
	var north = Vector2(0, 1)
	var west = Vector2(-1, 0)
	var south = Vector2(0, -1)
	var east = Vector2(1, 0)
	var ne = Vector2(1,-1)
	var nw = Vector2(-1,-1)
	var sw = Vector2(-1,1)
	var se = Vector2(1,1)
	
	assert Comparisons.is_counter_clockwise(north, west, origin)
	

		
		
func test_sort_events_by_angle():
	var a = VisibilityEvent.new(EVENT_TYPE.start_vertex, LineSegment.new(Vector2(1,0), Vector2(0,0)))
	var b = VisibilityEvent.new(EVENT_TYPE.end_vertex, LineSegment.new(Vector2(1,0),Vector2(0,0)))
	assert sort_events_by_angle(Vector2(0,0), [a,b]) == [b,a]

	var c = VisibilityEvent.new(EVENT_TYPE.end_vertex, LineSegment.new(Vector2(1,1),Vector2(0,0)))
#	assert sort_events_by_angle(Vector2(0,0), [a,b,c]) == [c,b,a]
#	assert sort_events_by_angle(Vector2(0,0), [a,c,b]) == [c,b,a]
#	assert sort_events_by_angle(Vector2(0,0), [b,a,c]) == [c,b,a]
#	assert sort_events_by_angle(Vector2(0,0), [b,c,a]) == [c,b,a]
#	assert sort_events_by_angle(Vector2(0,0), [c,a,b]) == [c,b,a]
#	assert sort_events_by_angle(Vector2(0,0), [c,b,a]) == [c,b,a]
	var dim = 250
	var top_left = Vector2(-dim, -dim)
	var top_right = Vector2(dim, -dim)
	var bottom_right = Vector2(dim,dim)
	var bottom_left = Vector2(-dim, dim)
	var origin = Vector2(0,0)
	
	
	var left = LineSegment.new(bottom_left, top_left)
	var bottom = LineSegment.new(bottom_right, bottom_left)
	var right =	LineSegment.new(top_right, bottom_right)
	var top = LineSegment.new(top_left, top_right)
	
	var revtop = reverse_segment(top)
	var revleft = reverse_segment(left)
	var revright = reverse_segment(right)
	var revbottom = reverse_segment(bottom)
	var line_left = Vector2(-40,50)
	var line_right = Vector2(40,50)


	var line = LineSegment.new(line_right, line_left)
	var rev_line = LineSegment.new(line_left, line_right)

	var names = {
		left: 'left',
		revleft: 'left_reverse',
		top: 'top',
		revtop: 'top_reverse',
		bottom: 'bottom',
		revbottom: 'bottom_reverse',
		right: 'right',
		revright: 'right_reverse',
		line: 'line',
		rev_line : 'line_reverse',
	}
	
	var left_start = VisibilityEvent.new(EVENT_TYPE.start_vertex, left)
	var left_end = 	VisibilityEvent.new(EVENT_TYPE.end_vertex, reverse_segment(left))
	var bot_start =	VisibilityEvent.new(EVENT_TYPE.start_vertex, bottom)
	var bot_end = VisibilityEvent.new(EVENT_TYPE.end_vertex, reverse_segment(bottom))
	var top_start = VisibilityEvent.new(EVENT_TYPE.start_vertex, top)
	var top_end = VisibilityEvent.new(EVENT_TYPE.end_vertex, reverse_segment(top))
	var right_start = VisibilityEvent.new(EVENT_TYPE.start_vertex, right)
	var right_end = VisibilityEvent.new(EVENT_TYPE.end_vertex, reverse_segment(right))
	var line_start = VisibilityEvent.new(EVENT_TYPE.start_vertex, line)
	var line_end = VisibilityEvent.new(EVENT_TYPE.end_vertex, rev_line)
	
	var expected_sort_order = [
		line_end, line_end,
		bot_end, left_start,
		left_end, top_start, 
		top_end, right_start, 
		right_end, bot_start,
		line_start, line_start,
	]
		
	
	var events = [left_start, left_end,
			bot_start, bot_end,
			top_start, top_end,
			right_start, right_end,
			line_start, line_end,
			line_start, line_end,
	]
	
	test_angle_sort_on_collection(origin, events, expected_sort_order)
	
	events = [
		right_start,
		line_start,
		line_start, 
		line_end,
		left_end,
		line_end,
		top_start, 
		bot_start,
		right_end,
		top_end,
		left_start, 
		bot_end,
	]
	
	test_angle_sort_on_collection(origin, events, expected_sort_order)


	
	
func test_angle_sort_on_collection(origin, collection, expected):
	var sorted = sort_events_by_angle(origin, collection)
	var act_segs = []
	var exp_segs = []
	for i in range(len(sorted)):
		exp_segs.append(expected[i].type)
		exp_segs.append(expected[i].segment.a)
		exp_segs.append(expected[i].segment.b)
		#exp_segs.append(names[expected[i].segment])
		act_segs.append(sorted[i].type)
		act_segs.append(sorted[i].segment.a)
		act_segs.append(sorted[i].segment.b)
		#act_segs.append(names[sorted[i].segment])
	
	print('done populating editor friendly data')
	for i in range(len(exp_segs)):
		assert exp_segs[i] == act_segs[i]
	
	
func test_storing_segments_as_events():
	pass


func test_calculate_visibility_polygon_with_no_line_segments():
	var segments = []
	var polygon = get_fov_from_segments(Vector2(0,0), segments)
	assert len(polygon) == 0

func test_store_segments_as_events():
	var dim = 250
	var top_left = Vector2(-dim, -dim)
	var top_right = Vector2(dim, -dim)
	var bottom_right = Vector2(dim,dim)
	var bottom_left = Vector2(-dim, dim)
	var origin = Vector2(0,0)
	
	
	var left = LineSegment.new(bottom_left, top_left)
	var bottom = LineSegment.new(bottom_right, bottom_left)
	var right =	LineSegment.new(top_right, bottom_right)
	var top = LineSegment.new(top_left, top_right)
	
	var revtop = reverse_segment(top)
	var revleft = reverse_segment(left)
	var revright = reverse_segment(right)
	var revbottom = reverse_segment(bottom)
	var segments = [
		top,
		bottom,
		right,
		left,
	]
	
	
	var data = store_segments_as_events(segments, Vector2(0,0))
	assert 'events' in data.keys()
	assert 'state' in data.keys()
	assert len(data['state'].data) == 1
	var state = data['state']
	var expected_seed_val = bottom
	var actual_seed = state.first()
	assert actual_seed.a == expected_seed_val.a
	assert actual_seed.b == expected_seed_val.b
	var events = data['events']
	

	
	var expected_events = [
		VisibilityEvent.new(EVENT_TYPE.start_vertex, top),
		VisibilityEvent.new(EVENT_TYPE.end_vertex, revtop),
		VisibilityEvent.new(EVENT_TYPE.start_vertex, bottom),
		VisibilityEvent.new(EVENT_TYPE.end_vertex, revbottom),
		VisibilityEvent.new(EVENT_TYPE.start_vertex, right),
		VisibilityEvent.new(EVENT_TYPE.end_vertex, revright),
		VisibilityEvent.new(EVENT_TYPE.start_vertex, left),	
		VisibilityEvent.new(EVENT_TYPE.end_vertex, revleft),
	]
	assert len(expected_events) == len(events)
	
	for idx in range(len(events)):
		var expect = expected_events[idx]
		var actual = events[idx]
		assert expect.type == actual.type
		assert expect.segment.a == actual.segment.a
		assert expect.segment.b == actual.segment.b
		
		

	var line_left = Vector2(-50,50)
	var line_right = Vector2(50,50)

	var line = LineSegment.new(line_right, line_left)
	var rev_line = LineSegment.new(line_left, line_right)
	
	
	events = store_segments_as_events( [left, bottom, top, right, line, rev_line], origin)['events']
	
	expected_events = [
		VisibilityEvent.new(EVENT_TYPE.start_vertex, left),	
		VisibilityEvent.new(EVENT_TYPE.end_vertex, reverse_segment(left)),
		VisibilityEvent.new(EVENT_TYPE.start_vertex, bottom),
		VisibilityEvent.new(EVENT_TYPE.end_vertex, reverse_segment(bottom)),
		VisibilityEvent.new(EVENT_TYPE.start_vertex, top),
		VisibilityEvent.new(EVENT_TYPE.end_vertex, reverse_segment(top)),
		VisibilityEvent.new(EVENT_TYPE.start_vertex, right),
		VisibilityEvent.new(EVENT_TYPE.end_vertex, reverse_segment(right)),
		VisibilityEvent.new(EVENT_TYPE.start_vertex, line),
		VisibilityEvent.new(EVENT_TYPE.end_vertex, rev_line),
		VisibilityEvent.new(EVENT_TYPE.start_vertex, line),
		VisibilityEvent.new(EVENT_TYPE.end_vertex, rev_line),
	]
	var ex_segs = []
	var act_segs = []
	for i in range(len(events)):
		ex_segs.append(expected_events[i].segment.a)
		ex_segs.append(expected_events[i].segment.b)
		act_segs.append(events[i].segment.a)
		act_segs.append(events[i].segment.b)
		
	assert len(events) == len(expected_events)
	for i in range(len(expected_events)):
		var actual = events[i]
		var expected = expected_events[i]
		
		assert actual.segment.a == expected.segment.a
		assert actual.segment.b == expected.segment.b
		assert actual.type == expected.type
	
	


	
func test_calculate_visibility_polygon_with_no_obstacle_only_boundary():
	
#	*-------------------------------*
#	|								|
#	|								|
#	|								|
#	|								|
#	|								|
#	|				X				|
#	|								|
#	|								|
#	|								|
#	|								|
#	|								|
#	*-------------------------------*
	var dim = 250
	var top_left = Vector2(-dim, -dim)
	var top_right = Vector2(dim, -dim)
	var bottom_right = Vector2(dim,dim)
	var bottom_left = Vector2(-dim, dim)
	var origin = Vector2(0,0)
	
	var left = LineSegment.new(bottom_left, top_left)
	var bottom = LineSegment.new(bottom_right, bottom_left)
	var right =	LineSegment.new(top_right, bottom_right)
	var top = LineSegment.new(top_left, top_right)
	var data = store_segments_as_events([top, right, left, bottom], origin)
	
	var events = sort_events_by_angle(origin, data['events'])

	var poly = find_visibility_polygon(events, data['state'], origin)
	assert len(poly) == 8
	assert Comparisons.new(origin).approx_equal_vector(poly[0], bottom_left)
	assert Comparisons.new(origin).approx_equal_vector(poly[1], bottom_left)
	assert Comparisons.new(origin).approx_equal_vector(poly[2], top_left)
	assert Comparisons.new(origin).approx_equal_vector(poly[3], top_left)
	assert Comparisons.new(origin).approx_equal_vector(poly[4], top_right)
	assert Comparisons.new(origin).approx_equal_vector(poly[5], top_right)
	assert Comparisons.new(origin).approx_equal_vector(poly[6], bottom_right)
	assert Comparisons.new(origin).approx_equal_vector(poly[7], bottom_right)

	
	
func test_calculate_visibility_polygon_with_a_polyline_obstacle():
	
#	*-------------------------------*
#	|								|
#	|								|
#	|								|
#	|								|
#	|								|
#	|				X				|
#	|								|
#	|			*-------*			|
#	|								|
#	|								|
#	|								|
#	*-------------------------------*
	# points
	var dim = 250
	var top_left = Vector2(-dim, -dim)
	var top_right = Vector2(dim, -dim)
	var bottom_right = Vector2(dim,dim)
	var bottom_left = Vector2(-dim, dim)
	var origin = Vector2(0,0)
	var line_left = Vector2(-50,50)
	var line_right = Vector2(50,50)

	#linesegments
	var left = LineSegment.new(bottom_left, top_left)
	var bottom = LineSegment.new(bottom_right, bottom_left)
	var right =	LineSegment.new(top_right, bottom_right)
	var top = LineSegment.new(top_left, top_right)
	var line = LineSegment.new(line_left, line_right)
	var rev_line = LineSegment.new(line_right, line_left)

	var segments = [top, right, left, bottom, line, rev_line]
	var data = store_segments_as_events(segments, origin)
	var events = sort_events_by_angle(origin, data['events'])
	var state = data['state']
	var poly = find_visibility_polygon(events, state, origin)
	assert len(poly) == 10
	assert Comparisons.new(origin).approx_equal_vector(poly[0], line_left)
	assert Comparisons.new(origin).approx_equal_vector(poly[1], bottom_left)
	assert Comparisons.new(origin).approx_equal_vector(poly[2], top_left)
	assert Comparisons.new(origin).approx_equal_vector(poly[3], top_left)
	assert Comparisons.new(origin).approx_equal_vector(poly[4], top_right)
	assert Comparisons.new(origin).approx_equal_vector(poly[5], top_right)
	assert Comparisons.new(origin).approx_equal_vector(poly[6], bottom_right)
	assert Comparisons.new(origin).approx_equal_vector(poly[7], bottom_right)
	assert Comparisons.new(origin).approx_equal_vector(poly[8], bottom_right)
	assert Comparisons.new(origin).approx_equal_vector(poly[9], line_right)


	
	var clean_poly = get_fov_from_segments(segments, origin)
	assert len(clean_poly) == 6
	assert Comparisons.new(origin).approx_equal_vector(clean_poly[0], line_left)
	assert Comparisons.new(origin).approx_equal_vector(clean_poly[1], bottom_left)
	assert Comparisons.new(origin).approx_equal_vector(clean_poly[2], top_left)
	assert Comparisons.new(origin).approx_equal_vector(clean_poly[3], top_right)
	assert Comparisons.new(origin).approx_equal_vector(clean_poly[4], bottom_right)
	assert Comparisons.new(origin).approx_equal_vector(clean_poly[5], line_right)
	
func test_calculate_visibility_polygon_with_a_convex_polygon_obstacle():
	
#	*-------------------------------*
#	|								|
#	|								|
#	|								|
#	|								|
#	|								|
#	|				X				|
#	|								|
#	|			 *-----*			|
#	|			 |	   |			|
#	|			 *-----*			|
#	|								|
#	*-------------------------------*
	var dim = 5
	var top_left = Vector2(-dim, -dim)
	var top_right = Vector2(dim, -dim)
	var bottom_right = Vector2(dim,dim)
	var bottom_left = Vector2(-dim, dim)
	var origin = Vector2(0,0)
	var box_dim = 2
	var box_top_left = Vector2(-box_dim,box_dim)
	var box_top_right = Vector2(box_dim,box_dim)
	var box_bottom_left = Vector2(-box_dim,box_dim+1)
	var box_bottom_right = Vector2(box_dim,box_dim+1)
	
	#linesegments
	var left = LineSegment.new(bottom_left, top_left)
	var bottom = LineSegment.new(bottom_right, bottom_left)
	var right =	LineSegment.new(top_right, bottom_right)
	var top = LineSegment.new(top_left, top_right)
	var box_left = LineSegment.new(box_top_left, box_bottom_left)
	var box_bottom = LineSegment.new(box_bottom_left, box_bottom_right)
	var box_right = LineSegment.new(box_bottom_right, box_top_right)
	var box_top = LineSegment.new(box_top_right, box_top_left)
	
	var segments = [
		box_right,   
		box_left,
		box_bottom, 
		box_top,
		top,
		right,
		left, 
		bottom, 
	] 
	var clean_poly = get_fov_from_segments(segments, origin)
	assert len(clean_poly) == 6
	assert Comparisons.new(origin).approx_equal_vector(clean_poly[0], box_top_left)
	assert Comparisons.new(origin).approx_equal_vector(clean_poly[1], bottom_left)
	assert Comparisons.new(origin).approx_equal_vector(clean_poly[2], top_left)
	assert Comparisons.new(origin).approx_equal_vector(clean_poly[3], top_right)
	assert Comparisons.new(origin).approx_equal_vector(clean_poly[4], bottom_right)
	assert Comparisons.new(origin).approx_equal_vector(clean_poly[5], box_top_right)
	
	test_near_block_at_angle()
	
	
func test_near_block_at_angle():
#	*-------------------------------*
#	|								|
#	|								|
#	|								|
#	|								|
#	|								|
#	|				X				|
#	|			   *-*-*			|
#	|			   *-*-*			|
#	|			 					|
#	|			 					|
#	|								|
#	*-------------------------------*
	var dim = 250.0
	var top_left = Vector2(-dim, -dim)
	var top_right = Vector2(dim, -dim)
	var bottom_right = Vector2(dim,dim)
	var bottom_left = Vector2(-dim, dim)
	var origin = Vector2(0.0,0.0)
	var box_top_left = Vector2(-.1,0.1)
	var box_top_right = Vector2(0.1,0.1)
	var box_bottom_left = Vector2(-0.1,0.2)
	var box_bottom_right = Vector2(0.1,0.2)
	var box_2_top_left = Vector2(-0.2, 0.1)
	var box_2_top_right = box_top_left
	var box_2_bottom_right = box_bottom_left
	var box_2_bottom_left = Vector2(-0.2, 0.2)
	
	#linesegments
	var left = LineSegment.new(bottom_left, top_left)
	var bottom = LineSegment.new(bottom_right, bottom_left)
	var right =	LineSegment.new(top_right, bottom_right)
	var top = LineSegment.new(top_left, top_right)
	var box_left = LineSegment.new(box_top_left, box_bottom_left)
	var box_bottom = LineSegment.new(box_bottom_left, box_bottom_right)
	var box_right = LineSegment.new(box_bottom_right, box_top_right)
	var box_top = LineSegment.new(box_top_right, box_top_left)
	var box_2_left = LineSegment.new(box_2_top_left, box_2_bottom_left)
	var box_2_bottom = LineSegment.new(box_2_bottom_left, box_2_bottom_right)
	var box_2_right = LineSegment.new(box_2_bottom_right, box_2_top_right)
	var box_2_top = LineSegment.new(box_2_top_right, box_2_top_left)
	
	var segments = [top, right, left, bottom, box_top, box_right, box_bottom, box_left, box_2_left, box_2_bottom, box_2_top, box_2_right, ]
	var clean_poly = get_fov_from_segments(segments, origin)
	assert Comparisons.new(origin).approx_equal_vector(clean_poly[0], box_2_top_right)
	assert Comparisons.new(origin).approx_equal_vector(clean_poly[1], box_2_top_left)
	assert Comparisons.new(origin).approx_equal_vector(clean_poly[2], Vector2(-dim, dim/2))
	assert Comparisons.new(origin).approx_equal_vector(clean_poly[3], top_left)
	assert Comparisons.new(origin).approx_equal_vector(clean_poly[4], top_right)
	assert Comparisons.new(origin).approx_equal_vector(clean_poly[5], bottom_right)
	assert Comparisons.new(origin).approx_equal_vector(clean_poly[6], box_top_right)
	
#
func test_calculate_visibility_polygon_with_a_concave_polygon_obstacle():
#	*-------------------------------*
#	|								|
#	|								|
#	|								|
#	|								|
#	|								|
#	|				X				|
#	|			*		*			|
#	|			\ \   / /			|
#	|			 \	*  / 			|
#	|			  \   /				|
#	|				*				|
#	*-------------------------------*
	
	var dim = 250
	var top_left = Vector2(-dim, -dim)
	var top_right = Vector2(dim, -dim)
	var bottom_right = Vector2(dim,dim)
	var bottom_left = Vector2(-dim, dim)
	var origin = Vector2(0,0)
	var box_top_left = Vector2(-50,50)
	var box_top_right = Vector2(50,50)
	var box_upper_mid = Vector2(0,100)
	var box_lower_mid = Vector2(0,200)
	
	#linesegments
	var left = LineSegment.new(bottom_left, top_left)
	var bottom = LineSegment.new(bottom_right, bottom_left)
	var right =	LineSegment.new(top_right, bottom_right)
	var top = LineSegment.new(top_left, top_right)
	var box_left_up = LineSegment.new(box_top_left, box_upper_mid)
	var box_left_down = LineSegment.new(box_top_left, box_lower_mid)
	var box_right_down = LineSegment.new(box_upper_mid, box_top_right)
	var box_right_up = LineSegment.new(box_lower_mid, box_top_right)
	
	var segments = [top, right, left, bottom, box_left_up, box_left_down, box_right_down, box_right_up]
	var clean_poly = get_fov_from_segments(segments, origin)
	assert len(clean_poly) == 7
	assert Comparisons.new(origin).approx_equal_vector(clean_poly[0], box_upper_mid)
	assert Comparisons.new(origin).approx_equal_vector(clean_poly[1], box_top_left)
	assert Comparisons.new(origin).approx_equal_vector(clean_poly[2], bottom_left)
	assert Comparisons.new(origin).approx_equal_vector(clean_poly[3], top_left)
	assert Comparisons.new(origin).approx_equal_vector(clean_poly[4], top_right)
	assert Comparisons.new(origin).approx_equal_vector(clean_poly[5], bottom_right)
	assert Comparisons.new(origin).approx_equal_vector(clean_poly[6], box_top_right)
	
