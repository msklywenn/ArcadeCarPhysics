using System;
using UnityEngine;

// TODO remove all reference to masses, work with accelerations directly

public class ArcadeCar : MonoBehaviour
{
    public LayerMask CollisionLayers;

    const float wheelWidth = 0.085f;

    public struct WheelData
    {
        // is wheel touched ground or not ?
        public bool isOnGround;

        // wheel ground touch point
        public RaycastHit touchPoint;

        // real yaw, after Ackermann steering correction
        public float yawRad;

        // visual rotation
        public float visualRotationRad;

        // suspension compression
        public float compression;

        // suspension compression on previous update
        public float compressionPrev;

        public float rotationsPerSecond;
        public float rotationsPerSecondSpeed;
    }

    [Serializable]
    public class Axle
    {
        [Header("Axle settings")]
        [Tooltip("Axle width")]
        public float width = 0.4f;

        [Tooltip("Axle offset")]
        public Vector2 offset = Vector2.zero;

        [Tooltip("Current steering angle (in degrees)")]
        public float steerAngle = 0.0f;

        [Tooltip("Is axle powered by engine")]
        public bool isPowered = false;

        //--
        [Header("Wheel settings")]
        [Tooltip("Wheel radius in meters")]
        public float radius = 0.3f;

        [Range(0.0f, 1.0f)]
        [Tooltip("Tire laterial friction normalized to 0..1")]
        public float laterialFriction = 0.1f;

        [Range(0.0f, 0.1f)]
        [Tooltip("Rolling friction, normalized to 0..1")]
        public float rollingFriction = 0.01f;

        [Tooltip("Brake force magnitude")]
        public float brakeForceMag = 4.0f;

        //--
        [Header("Suspension settings")]
        [Tooltip("Suspension Stiffness (Suspension 'Power'")]
        public float stiffness = 8500.0f;

        [Tooltip("Suspension Damping (Suspension 'Bounce')")]
        public float damping = 3000.0f;

        [HideInInspector]
        [Tooltip("Suspension Restitution (Not used now)")]
        public float restitution = 1.0f;

        [Tooltip("Relaxed suspension length")]
        public float lengthRelaxed = 0.55f;

        [Tooltip("Stabilizer bar anti-roll force")]
        public float antiRollForce = 10000.0f;

        [HideInInspector]
        public WheelData wheelDataL;

        [HideInInspector]
        public WheelData wheelDataR;

        //--
        [Header("Visual settings")]
        [Tooltip("Visual scale for wheels")]
        public float visualScale = 0.03270531f;

        [Tooltip("Wheel actor left")]
        public GameObject wheelVisualLeft;

        [Tooltip("Wheel actor right")]
        public GameObject wheelVisualRight;

        //--
        [Header("Debug settings")]
        [Tooltip("Debug color for axle")]
        public Color debugColor = Color.white;
    }

    public Vector3 centerOfMass = Vector3.zero;

    [Tooltip("Only used as a reference for downforce and steering settings")]
    public float MaxSpeed = 200f; // used only for settings

    [Header("Engine")]
    public float ForwardTopSpeed = 130f;
    [Tooltip("How quickly to reach top speed. >1 is slowly, <1 is quickly")]
    [Min(0.01f)] public float ForwardAccelerationFactor = 1f / 3f;
    [Tooltip("Time to reach top speed")]
    [Min(0.01f)] public float ForwardAccelerationTime = 10f;

    public float ReverseTopSpeed = 40f;
    [Min(0.01f)] public float ReverseAccelerationFactor = 1f / 3f;
    [Min(0.01f)] public float ReverseAccelerationTime = 10f;

    [Header("Steering")]
    public float steeringSpeed = 60f;
    public float steeringResetSpeed = 50f;
    // x - speed in km/h
    // y - angle in degrees
    public float SlowSteerAngleLimit = 40f;
    public float FastSteerAngleLimit = 4f;
    public float SlowFastSteerCurve = 0.9f;

    [Header("Other")]
    [Tooltip("Stabilization in flight (torque)")]
    public float flightStabilizationForce = 8.0f;

    [Tooltip("Stabilization in flight (Ang velocity damping)")]
    public float flightStabilizationDamping = 0.0f;

    public float DownForceIntensity = 0.7f;

    [Tooltip("Downforce")]
    public float DownForce = 10f;

    [Tooltip("Rotation speed when braking and accelerating at the same time")]
    public float burnRotationSpeed = 5f;
    public float wheelRotationSmoothing = 0.1f;

    [Header("Axles")]
    public Axle[] axles = new Axle[2];

    [Header("Input")]
    public bool controllable = true;
    [Range(0, 1)] public float SteeringDeadzone = 0.01f;

    [Header("Debug")]
    public bool debugDraw = true;

    //////////////////////////////////////////////////////////////////////////////////////////////////////

    bool isBrake = false;
    bool isHandBrake = false;
    bool isAcceleration = false;
    bool isReverseAcceleration = false;
    float accelerationForceMagnitude = 0.0f;
    Rigidbody rb = null;

    // UI style for debug render
    static GUIStyle style = new GUIStyle();

    // For alloc-free raycasts
    Ray wheelRay = new Ray();
    RaycastHit[] wheelRayHits = new RaycastHit[4];

    void Teleport(Vector3 position)
    {
        position += new Vector3(UnityEngine.Random.Range(-1.0f, 1.0f), 0.0f, UnityEngine.Random.Range(-1.0f, 1.0f));
        float yaw = transform.eulerAngles.y + UnityEngine.Random.Range(-10.0f, 10.0f);

        transform.SetPositionAndRotation(position, Quaternion.Euler(new Vector3(0.0f, yaw, 0.0f)));

        rb.velocity = new Vector3(0f, 0f, 0f);
        rb.angularVelocity = new Vector3(0f, 0f, 0f);

        for (int axleIndex = 0; axleIndex < axles.Length; axleIndex++)
            axles[axleIndex].steerAngle = 0.0f;
    }

    void Start()
    {
        style.normal.textColor = Color.red;
        
        TryGetComponent(out rb);
        rb.centerOfMass = centerOfMass;
    }

    void OnValidate()
    {
        //HACK: to apply steering in editor
        if (rb == null)
            TryGetComponent(out rb);

        ApplyVisual();
        CalculateAckermannSteering();
    }

    static float EasyCurve(float t, float intensity)
    {
        return Mathf.Pow(MathF.Sin(Mathf.Clamp01(t) * MathF.PI / 2f), intensity);
    }

    static float EasyCurveInverse(float t, float intensity)
    {
        return MathF.Asin(MathF.Pow(Mathf.Clamp01(t), 1f / intensity)) * 2f / MathF.PI;
    }

    float GetAccelerationForceMagnitude(float topSpeed, float factor, float duration, float speedMetersPerSec, float dt)
    {
        float speedKmH = speedMetersPerSec * 3.6f;

        float mass = rb.mass;

        float t = duration * EasyCurveInverse(speedKmH / topSpeed, factor);
        float desiredSpeed = topSpeed * EasyCurve((t + dt) / duration, factor);
        desiredSpeed = Mathf.Clamp(desiredSpeed, 0f, topSpeed);
        float acceleration = desiredSpeed - speedKmH;
        
        acceleration /= 3.6f; //to meters per sec
        float intensity = acceleration * mass;
        intensity = Mathf.Max(intensity, 0.0f);
        return intensity;
    }

    public float GetSpeed()
    {
        Vector3 velocity = rb.velocity;

        Vector3 wsForward = rb.transform.rotation * Vector3.forward;
        float vProj = Vector3.Dot(velocity, wsForward);
        Vector3 projVelocity = vProj * wsForward;
        float speed = projVelocity.magnitude * Mathf.Sign(vProj);
        return speed;
    }

    float CalcAccelerationForceMagnitude()
    {
        if (!isAcceleration && !isReverseAcceleration)
            return 0.0f;

        float speed = GetSpeed();
        float dt = Time.fixedDeltaTime;

        if (isAcceleration)
            return GetAccelerationForceMagnitude(ForwardTopSpeed, ForwardAccelerationFactor, ForwardAccelerationTime, speed, dt);
        else
            return -GetAccelerationForceMagnitude(ReverseTopSpeed, ReverseAccelerationFactor, ReverseAccelerationTime, -speed, dt);
    }

    void Steering(float steeringWheel, float speed)
    {
        float speedKph = Mathf.Abs(speed) * 3.6f;
        float steering;
        if (Mathf.Abs(steeringWheel) > SteeringDeadzone)
        {
            steering = steeringWheel * steeringSpeed * Time.fixedDeltaTime;
        }
        else
        {
            float resetSign = -Mathf.Sign(axles[0].steerAngle);
            float reset = Mathf.Min(Mathf.Abs(axles[0].steerAngle), steeringResetSpeed * Time.fixedDeltaTime);
            steering = reset * resetSign;
        }

        float newSteerAngle = axles[0].steerAngle + steering;
        float sgn = Mathf.Sign(newSteerAngle);
        float steerLimit = Mathf.Lerp(SlowSteerAngleLimit, FastSteerAngleLimit, 
            EasyCurve(speedKph / MaxSpeed, SlowFastSteerCurve));
        newSteerAngle = Mathf.Min(Math.Abs(newSteerAngle), steerLimit) * sgn;
        axles[0].steerAngle = newSteerAngle;
    }

    void UpdateInput()
    {
        float v = 0f;
        float wheel = 0f;
        bool isBrakeNow = false;
        bool isHandBrakeNow = false;

        if (controllable)
        {
            v = Input.GetAxis("Vertical");
            wheel = Input.GetAxis("Horizontal");
            if (Input.GetKey(KeyCode.R)) 
                ResetToValidPosition();
            isHandBrakeNow = Input.GetKey(KeyCode.Space);
            isBrakeNow = Input.GetKey(KeyCode.RightControl) || Input.GetKey(KeyCode.LeftControl);
        }

        float speed = GetSpeed();
        isAcceleration = v > 0.4f && speed > -0.5f;
        isReverseAcceleration = v < -0.4f && speed < 0.5f;

        isBrake = isBrakeNow;

        // hand brake + acceleration = power slide
        isHandBrake = isHandBrakeNow && !isAcceleration && !isReverseAcceleration;

        if (controllable)
            Steering(wheel, speed);
    }

    private void ResetToValidPosition()
    {
        Debug.Log("Reset pressed");
        Ray resetRay = new Ray();

        // trace top-down
        resetRay.origin = transform.position + new Vector3(0.0f, 100.0f, 0.0f);
        resetRay.direction = new Vector3(0.0f, -1.0f, 0.0f);

        RaycastHit[] resetRayHits = new RaycastHit[16];

        int numHits = Physics.RaycastNonAlloc(resetRay, resetRayHits, 250.0f);

        if (numHits > 0)
        {
            float nearestDistance = float.MaxValue;
            for (int j = 0; j < numHits; j++)
            {
                if (resetRayHits[j].collider != null && resetRayHits[j].collider.isTrigger)
                {
                    // skip triggers
                    continue;
                }

                // Filter contacts with car body
                if (resetRayHits[j].rigidbody == rb)
                {
                    continue;
                }

                if (resetRayHits[j].distance < nearestDistance)
                {
                    nearestDistance = resetRayHits[j].distance;
                }
            }

            // -4 meters from surface
            nearestDistance -= 4.0f;
            Vector3 resetPos = resetRay.origin + resetRay.direction * nearestDistance;
            Teleport(resetPos);
        }
        else
        {
            // Hard reset
            Teleport(new Vector3(-69.48f, 5.25f, 132.71f));
        }
    }

    void Update()
    {
        ApplyVisual();
    }

    void FixedUpdate()
    {
        UpdateInput();

        accelerationForceMagnitude = CalcAccelerationForceMagnitude();

        CalculateAckermannSteering();

        int numberOfPoweredWheels = 0;
        for (int axleIndex = 0; axleIndex < axles.Length; axleIndex++)
            if (axles[axleIndex].isPowered)
                numberOfPoweredWheels += 2;

        int totalWheelsCount = axles.Length * 2;
        for (int axleIndex = 0; axleIndex < axles.Length; axleIndex++)
            CalculateAxleForces(axles[axleIndex], totalWheelsCount, numberOfPoweredWheels);

        bool inAir = AreAllWheelsInAir();

        if (inAir)
            KeepUpwards();
        else
            Downforce();
    }

    private bool AreAllWheelsInAir()
    {
        bool allWheelIsOnAir = true;
        for (int axleIndex = 0; axleIndex < axles.Length; axleIndex++)
        {
            if (axles[axleIndex].wheelDataL.isOnGround || axles[axleIndex].wheelDataR.isOnGround)
            {
                allWheelIsOnAir = false;
                break;
            }
        }

        return allWheelIsOnAir;
    }

    private void KeepUpwards()
    {
        // Try to keep vehicle parallel to the ground while jumping
        Vector3 carUp = transform.TransformDirection(new Vector3(0.0f, 1.0f, 0.0f));
        Vector3 worldUp = new Vector3(0.0f, 1.0f, 0.0f);

        // Flight stabilization from
        // TODO: this is probably GPL code so replace with my previous solution
        // https://github.com/supertuxkart/stk-code/blob/master/src/physics/btKart.cpp#L455

        // Length of axis depends on the angle - i.e. the further awat
        // the kart is from being upright, the larger the applied impulse
        // will be, resulting in fast changes when the kart is on its
        // side, but not overcompensating (and therefore shaking) when
        // the kart is not much away from being upright.
        Vector3 axis = Vector3.Cross(carUp, worldUp);
        //axis.Normalize ();

        float mass = rb.mass;

        // angular velocity damping
        Vector3 angVel = rb.angularVelocity;

        Vector3 angVelDamping = angVel;
        angVelDamping.y = 0.0f;
        angVelDamping *= Mathf.Clamp01(flightStabilizationDamping * Time.fixedDeltaTime);

        //Debug.Log(string.Format("Ang {0}, Damping {1}", angVel, angVelDamping));
        rb.angularVelocity = angVel - angVelDamping;

        // in flight roll stabilization
        rb.AddTorque(flightStabilizationForce * mass * axis);
    }

    private void Downforce()
    {
        // downforce
        Vector3 carDown = transform.TransformDirection(new Vector3(0.0f, -1.0f, 0.0f));

        float speed = GetSpeed();
        float speedKmH = Mathf.Abs(speed) * 3.6f;

        float downForceAmount = 1f - EasyCurve(1f - speedKmH / MaxSpeed, DownForceIntensity);

        float mass = rb.mass;

        rb.AddForce(DownForce * downForceAmount * mass * carDown);

        //Debug.Log(string.Format("{0} downforce", downForceAmount * downForce));
    }

    void OnGUI()
    {
        if (!controllable)
            return;

        float speed = GetSpeed();
        float speedKmH = speed * 3.6f;
        GUI.Label(new Rect(30.0f, 20.0f, 150, 130), string.Format("{0:F2} km/h", speedKmH), style);

        float yPos = 60.0f;
        for (int axleIndex = 0; axleIndex < axles.Length; axleIndex++)
        {
            GUI.Label(new Rect(30.0f, yPos, 1500, 130),
                $"Axle {axleIndex}, steering angle {axles[axleIndex].steerAngle:F2}, compression: {axles[axleIndex].wheelDataL.compression:F2} {axles[axleIndex].wheelDataR.compression:F2}", style);
            yPos += 18.0f;
        }

        Camera cam = Camera.current;
        if (cam == null)
        {
            return;
        }

        if (debugDraw)
        {
            foreach (Axle axle in axles)
            {
                Vector3 localL = new Vector3(axle.width * -0.5f, axle.offset.y, axle.offset.x);
                Vector3 localR = new Vector3(axle.width * 0.5f, axle.offset.y, axle.offset.x);

                Vector3 wsL = transform.TransformPoint(localL);
                Vector3 wsR = transform.TransformPoint(localR);

                Vector3 screenPos = cam.WorldToScreenPoint(wsL);
                GUI.Label(new Rect(screenPos.x, Screen.height - screenPos.y, 150, 130),
                    axle.wheelDataL.compression.ToString("F2"), style);
                screenPos = cam.WorldToScreenPoint(wsR);
                GUI.Label(new Rect(screenPos.x, Screen.height - screenPos.y, 150, 130),
                    axle.wheelDataR.compression.ToString("F2"), style);
            }
        }
    }

    void AddForceAtPosition(in Vector3 force, in Vector3 position)
    {
        rb.AddForceAtPosition(force, position);
        //Debug.DrawRay(position, force, Color.magenta);
    }

    bool RayCast(in Ray ray, float maxDistance, ref RaycastHit nearestHit)
    {
        int numHits = Physics.RaycastNonAlloc(ray, wheelRayHits, maxDistance, CollisionLayers, QueryTriggerInteraction.Ignore);
        if (numHits <= 0)
            return false;
        bool hasHit = false;
        for (int i = 0; i < numHits; i++)
        {
            if (wheelRayHits[i].normal.y > 0.6f && (!hasHit || nearestHit.distance > wheelRayHits[i].distance))
            {
                nearestHit = wheelRayHits[i];
                hasHit = true;
            }
        }
        return hasHit && nearestHit.distance <= maxDistance;
    }

    void CalculateWheelForces(Axle axle, in Vector3 wsDownDirection, ref WheelData wheelData, in Vector3 wsAttachPoint, int totalWheelsCount, int numberOfPoweredWheels)
    {
        float dt = Time.fixedDeltaTime;

        // Get wheel world space rotation and axes
        Quaternion localWheelRot = Quaternion.Euler(new Vector3(0.0f, wheelData.yawRad * Mathf.Rad2Deg, 0.0f));
        Quaternion wsWheelRot = transform.rotation * localWheelRot;

        // Wheel axle left direction
        Vector3 wsAxleLeft = wsWheelRot * Vector3.left;

        wheelData.isOnGround = false;
        wheelRay.direction = wsDownDirection;

        // Ray cast (better to use shape cast here, but Unity does not support shape casts)
        // TODO: replace with a single SphereCast
        float traceLen = axle.lengthRelaxed + axle.radius;

        wheelRay.origin = wsAttachPoint + wsAxleLeft * wheelWidth;
        RaycastHit s1 = new RaycastHit();
        bool b1 = RayCast(wheelRay, traceLen, ref s1);

        wheelRay.origin = wsAttachPoint - wsAxleLeft * wheelWidth;
        RaycastHit s2 = new RaycastHit();
        bool b2 = RayCast(wheelRay, traceLen, ref s2);

        wheelRay.origin = wsAttachPoint;
        bool isCollided = RayCast(wheelRay, traceLen, ref wheelData.touchPoint);

        // No wheel contant found
        if (!isCollided || !b1 || !b2)
        {
            // wheel do not touch the ground (relaxing spring)
            float relaxSpeed = 1.0f;
            wheelData.compressionPrev = wheelData.compression;
            wheelData.compression = Mathf.Clamp01(wheelData.compression - dt * relaxSpeed);
            return;
        }

        // Consider wheel radius
        float suspLenNow = wheelData.touchPoint.distance - axle.radius;

        Debug.AssertFormat(suspLenNow <= traceLen, "Sanity check failed.");

        wheelData.isOnGround = true;

        //
        // Calculate suspension force
        //
        // Spring force - want's go to position 0
        // Damping force - want's go to velocity 0
        //
        ////////////////////////////////////////////////////////////////////////////////////////////////////

        float suspForceMag = 0.0f;

        // Positive value means that the spring is compressed
        // Negative value means that the spring is elongated.

        wheelData.compression = 1.0f - Mathf.Clamp01(suspLenNow / axle.lengthRelaxed);

        // Hooke's law (springs)
        // F = -k x 

        // Spring force (try to reset compression from spring)
        float springForce = wheelData.compression * -axle.stiffness;
        suspForceMag += springForce;

        // Damping force (try to reset velocity to 0)
        float suspCompressionVelocity = (wheelData.compression - wheelData.compressionPrev) / dt;
        wheelData.compressionPrev = wheelData.compression;

        float damperForce = -suspCompressionVelocity * axle.damping;
        suspForceMag += damperForce;

        // Only consider component of force that is along the contact normal.
        float denom = Vector3.Dot(wheelData.touchPoint.normal, -wsDownDirection);
        suspForceMag *= denom;

        // Apply suspension force
        Vector3 suspForce = wsDownDirection * suspForceMag;
        AddForceAtPosition(suspForce, wheelData.touchPoint.point);
        if (debugDraw)
            Debug.DrawRay(wheelData.touchPoint.point, suspForce / rb.mass, Color.yellow);

        //
        // Calculate friction forces
        //
        ////////////////////////////////////////////////////////////////////////////////////////////////////
        /// 

        Vector3 wheelVelocity = rb.GetPointVelocity(wheelData.touchPoint.point);

        // Contact basis (can be different from wheel basis)
        Vector3 c_up = wheelData.touchPoint.normal;
        Vector3 c_fwd = Vector3.ProjectOnPlane(wsWheelRot * Vector3.forward, c_up);

        // Calculate sliding velocity (velocity without normal force)
        Vector3 slideVelocity = Vector3.ProjectOnPlane(wheelVelocity, c_up);

        // Calculate current sliding force
        Vector3 slidingForce = rb.mass / dt / totalWheelsCount * slideVelocity;

        if (debugDraw)
            Debug.DrawRay(wheelData.touchPoint.point, slideVelocity / rb.mass, Color.red);

        float lateralFriction = Mathf.Clamp01(axle.laterialFriction);

        // Simulate perfect static friction
        Vector3 frictionForce = slidingForce * -lateralFriction;

        // Remove friction along roll-direction of wheel 
        Vector3 longitudinalForce = Vector3.Dot(frictionForce, c_fwd) * c_fwd;

        // Apply braking force or rolling resistance force or nothing
        bool stop = Mathf.Abs(GetSpeed()) < 0.1f;
        if (stop || isBrake || isHandBrake)
        {
            float mag = longitudinalForce.magnitude;
            float brakeForce = mag > 0f ? Mathf.Clamp(axle.brakeForceMag * rb.mass, 0.0f, mag) / mag : 0f;

            if (isHandBrake)
            {
                // hand brake are not powerful enough ;)
                brakeForce *= 0.8f;
            }

            longitudinalForce *= 1f - brakeForce;
        }
        else if (!isAcceleration && !isReverseAcceleration)
        {
            // Apply rolling-friction (automatic slow-down) only if player don't press the accelerator
            float rollingK = 1.0f - Mathf.Clamp01(axle.rollingFriction);
            longitudinalForce *= rollingK;
        }

        if (debugDraw)
        {
            Debug.DrawRay(wheelData.touchPoint.point, frictionForce / rb.mass, Color.red);
            Debug.DrawRay(wheelData.touchPoint.point, (frictionForce - longitudinalForce) / rb.mass, Color.blue);
            Debug.DrawRay(wheelData.touchPoint.point, longitudinalForce / rb.mass, Color.white);
        }

        frictionForce -= longitudinalForce;
            
        // Apply resulting force
        AddForceAtPosition(frictionForce, wheelData.touchPoint.point);

        // Engine force
        if (axle.isPowered && Mathf.Abs(accelerationForceMagnitude) > 0.01f)
        {
            Vector3 accForcePoint = wheelData.touchPoint.point - (wsDownDirection * 0.2f);
            Vector3 engineForce = c_fwd * accelerationForceMagnitude / numberOfPoweredWheels / dt;
            AddForceAtPosition(engineForce, accForcePoint);

            if (debugDraw)
                Debug.DrawRay(accForcePoint, engineForce / rb.mass, Color.green);
        }
    }

    void CalculateAxleForces(Axle axle, int totalWheelsCount, int numberOfPoweredWheels)
    {
        Vector3 wsDownDirection = transform.TransformDirection(Vector3.down);
        wsDownDirection.Normalize();

        Vector3 localL = new Vector3(axle.width * -0.5f, axle.offset.y, axle.offset.x);
        Vector3 localR = new Vector3(axle.width * 0.5f, axle.offset.y, axle.offset.x);

        Vector3 wsL = transform.TransformPoint(localL);
        Vector3 wsR = transform.TransformPoint(localR);

        CalculateWheelForces(axle, wsDownDirection, ref axle.wheelDataL, wsL, totalWheelsCount, numberOfPoweredWheels);
        CalculateWheelForces(axle, wsDownDirection, ref axle.wheelDataR, wsR, totalWheelsCount, numberOfPoweredWheels);

        // http://projects.edy.es/trac/edy_vehicle-physics/wiki/TheStabilizerBars
        // Apply "stablizer bar" forces
        float travelL = 1.0f - Mathf.Clamp01(axle.wheelDataL.compression);
        float travelR = 1.0f - Mathf.Clamp01(axle.wheelDataR.compression);

        float antiRollForce = (travelL - travelR) * axle.antiRollForce;
        if (axle.wheelDataL.isOnGround)
        {
            AddForceAtPosition(wsDownDirection * antiRollForce, axle.wheelDataL.touchPoint.point);
            if (debugDraw)
                Debug.DrawRay(axle.wheelDataL.touchPoint.point, wsDownDirection * antiRollForce / rb.mass, Color.magenta);
        }

        if (axle.wheelDataR.isOnGround)
        {
            AddForceAtPosition(wsDownDirection * -antiRollForce, axle.wheelDataR.touchPoint.point);
            if (debugDraw)
                Debug.DrawRay(axle.wheelDataR.touchPoint.point, wsDownDirection * -antiRollForce / rb.mass, Color.magenta);
        }
    }

    void CalculateAckermannSteering()
    {
        // Copy desired steering
        for (int axleIndex = 0; axleIndex < axles.Length; axleIndex++)
        {
            float steerAngleRad = axles[axleIndex].steerAngle * Mathf.Deg2Rad;

            axles[axleIndex].wheelDataL.yawRad = steerAngleRad;
            axles[axleIndex].wheelDataR.yawRad = steerAngleRad;
        }

        if (axles.Length != 2)
        {
            Debug.LogWarning("Ackermann work only for 2 axle vehicles.");
            return;
        }

        Axle frontAxle = axles[0];
        Axle rearAxle = axles[1];

        if (Mathf.Abs(rearAxle.steerAngle) > 0.0001f)
        {
            Debug.LogWarning("Ackermann work only for vehicles with forward steering axle.");
            return;
        }

        // Calculate our chassis (remove scale)
        Vector3 axleDiff = transform.TransformPoint(new Vector3(0.0f, frontAxle.offset.y, frontAxle.offset.x)) - transform.TransformPoint(new Vector3(0.0f, rearAxle.offset.y, rearAxle.offset.x));
        float axleSeparation = axleDiff.magnitude;

        Vector3 wheelDiff = transform.TransformPoint(new Vector3(frontAxle.width * -0.5f, frontAxle.offset.y, frontAxle.offset.x)) - transform.TransformPoint(new Vector3(frontAxle.width * 0.5f, frontAxle.offset.y, frontAxle.offset.x));
        float wheelsSeparation = wheelDiff.magnitude;

        // Get turning circle radius for steering angle input
        float turningCircleRadius = axleSeparation / Mathf.Tan(frontAxle.steerAngle * Mathf.Deg2Rad);

        // Make front inside tire turn sharper and outside tire less sharp based on turning circle radius
        float steerAngleLeft = Mathf.Atan(axleSeparation / (turningCircleRadius + (wheelsSeparation / 2)));
        float steerAngleRight = Mathf.Atan(axleSeparation / (turningCircleRadius - (wheelsSeparation / 2)));

        frontAxle.wheelDataL.yawRad = steerAngleLeft;
        frontAxle.wheelDataR.yawRad = steerAngleRight;
    }

    void CalculateWheelVisualTransform(in Vector3 wsAttachPoint, in Vector3 wsDownDirection, Axle axle, in WheelData data, bool leftWheel, out Vector3 pos, out Quaternion rot)
    {
        float suspCurrentLen = Mathf.Clamp01(1.0f - data.compression) * axle.lengthRelaxed;

        pos = wsAttachPoint + wsDownDirection * suspCurrentLen;

        float additionalYaw;
        float additionalMul;
        if (leftWheel)
        {
            additionalYaw = 180.0f;
            additionalMul = -Mathf.Rad2Deg;
        }
        else
        {
            additionalYaw = 0.0f;
            additionalMul = Mathf.Rad2Deg;
        }

        Quaternion localWheelRot = Quaternion.Euler(
            data.visualRotationRad * additionalMul,
            additionalYaw + data.yawRad * Mathf.Rad2Deg,
            0.0f);
        rot = transform.rotation * localWheelRot;
    }

    void CalculateWheelRotationFromSpeed(Axle axle, ref WheelData data, in Vector3 wsPos)
    {
        if (rb == null)
        {
            data.visualRotationRad = 0.0f;
            return;
        }

        float rps;
        if (axle.isPowered && isBrake && (isAcceleration || isReverseAcceleration))
        {
            rps = isAcceleration ? burnRotationSpeed : -burnRotationSpeed;
        }
        else
        {

            Quaternion localWheelRot = Quaternion.Euler(new Vector3(0.0f, data.yawRad * Mathf.Rad2Deg, 0.0f));
            Quaternion wsWheelRot = transform.rotation * localWheelRot;

            Vector3 wsWheelForward = wsWheelRot * Vector3.forward;
            Vector3 velocityQueryPos = data.isOnGround ? data.touchPoint.point : wsPos;
            Vector3 wheelVelocity = rb.GetPointVelocity(velocityQueryPos);

            // Longitudinal speed (meters/sec)
            float tireLongSpeed = Vector3.Dot(wheelVelocity, wsWheelForward);

            // Circle length = 2 * PI * R
            float wheelLengthMeters = 2 * Mathf.PI * axle.radius;

            // Wheel "Revolutions per second";
            rps = tireLongSpeed / wheelLengthMeters;
        }
            
        data.rotationsPerSecond = Mathf.SmoothDamp(data.rotationsPerSecond, rps,
            ref data.rotationsPerSecondSpeed, wheelRotationSmoothing);

        float deltaRot = Mathf.PI * 2.0f * data.rotationsPerSecond * Time.deltaTime;

        data.visualRotationRad += deltaRot;
    }

    void ApplyVisual()
    {
        Vector3 wsDownDirection = transform.TransformDirection(Vector3.down);
        wsDownDirection.Normalize();

        for (int axleIndex = 0; axleIndex < axles.Length; axleIndex++)
        {
            Axle axle = axles[axleIndex];

            Vector3 localL = new Vector3(axle.width * -0.5f, axle.offset.y, axle.offset.x);
            Vector3 localR = new Vector3(axle.width * 0.5f, axle.offset.y, axle.offset.x);

            Vector3 wsL = transform.TransformPoint(localL);
            Vector3 wsR = transform.TransformPoint(localR);

            Vector3 wsPos;
            Quaternion wsRot;

            if (axle.wheelVisualLeft != null)
            {
                CalculateWheelVisualTransform(wsL, wsDownDirection, axle, axle.wheelDataL, true, out wsPos, out wsRot);
                axle.wheelVisualLeft.transform.SetPositionAndRotation(wsPos, wsRot);
                axle.wheelVisualLeft.transform.localScale
                    = new Vector3(axle.radius, axle.radius, axle.radius) * axle.visualScale;

                CalculateWheelRotationFromSpeed(axle, ref axle.wheelDataL, wsPos);
            }

            if (axle.wheelVisualRight != null)
            {
                CalculateWheelVisualTransform(wsR, wsDownDirection, axle, axle.wheelDataR, false, out wsPos, out wsRot);
                axle.wheelVisualRight.transform.SetPositionAndRotation(wsPos, wsRot);
                axle.wheelVisualRight.transform.localScale
                    = new Vector3(axle.radius, axle.radius, axle.radius) * axle.visualScale;

                CalculateWheelRotationFromSpeed(axle, ref axle.wheelDataR, wsPos);
            }
        }
    }

#if UNITY_EDITOR
    void OnDrawGizmosAxle(Vector3 wsDownDirection, Axle axle)
    {
        Vector3 localL = new Vector3(axle.width * -0.5f, axle.offset.y, axle.offset.x);
        Vector3 localR = new Vector3(axle.width * 0.5f, axle.offset.y, axle.offset.x);

        Vector3 wsL = transform.TransformPoint(localL);
        Vector3 wsR = transform.TransformPoint(localR);

        Gizmos.color = axle.debugColor;

        //draw axle
        Gizmos.DrawLine(wsL, wsR);

        //draw line to com
        Gizmos.DrawLine(transform.TransformPoint(new Vector3(0.0f, axle.offset.y, axle.offset.x)), transform.TransformPoint(centerOfMass));

        DrawWheelGizmo(ref wsDownDirection, axle, true, axle.wheelDataL, wsL);
        DrawWheelGizmo(ref wsDownDirection, axle, false, axle.wheelDataL, wsL);
    }

    private void DrawWheelGizmo(ref Vector3 wsDownDirection, Axle axle, bool leftWheel, in WheelData wheelData, in Vector3 wsFrom)
    {
        Gizmos.color = wheelData.isOnGround ? Color.yellow : axle.debugColor;
        UnityEditor.Handles.color = Gizmos.color;

        float suspCurrentLen = Mathf.Clamp01(1.0f - wheelData.compression) * axle.lengthRelaxed;

        Vector3 wsTo = wsFrom + wsDownDirection * suspCurrentLen;

        // Draw suspension
        Gizmos.DrawLine(wsFrom, wsTo);

        Quaternion localWheelRot = Quaternion.Euler(new Vector3(0.0f, wheelData.yawRad * Mathf.Rad2Deg, 0.0f));
        Quaternion wsWheelRot = transform.rotation * localWheelRot;

        Vector3 localAxle = leftWheel ? Vector3.left : Vector3.right;
        Vector3 wsAxle = wsWheelRot * localAxle;
        Vector3 wsForward = wsWheelRot * Vector3.forward;

        // Draw wheel axle
        Gizmos.DrawLine(wsTo, wsTo + wsAxle * 0.1f);
        Gizmos.DrawLine(wsTo, wsTo + wsForward * axle.radius);

        // Draw wheel
        UnityEditor.Handles.DrawWireDisc(wsTo, wsAxle, axle.radius);

        UnityEditor.Handles.DrawWireDisc(wsTo + wsAxle * wheelWidth, wsAxle, axle.radius);
        UnityEditor.Handles.DrawWireDisc(wsTo - wsAxle * wheelWidth, wsAxle, axle.radius);
    }

    void OnDrawGizmos()
    {
        Vector3 wsDownDirection = transform.TransformDirection(Vector3.down);
        wsDownDirection.Normalize();
        foreach (Axle axle in axles)
            OnDrawGizmosAxle(wsDownDirection, axle);
        Gizmos.DrawSphere(transform.TransformPoint(centerOfMass), 0.1f);
    }
#endif
}

