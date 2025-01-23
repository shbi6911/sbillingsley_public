function vector = TransformFromInertialToBody(vector_inertial, euler_angles)

vector = rotation321(euler_angles)*vector_inertial;