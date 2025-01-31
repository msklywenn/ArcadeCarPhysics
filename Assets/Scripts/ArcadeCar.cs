using System;
using UnityEngine;

static class Utilities
{
    public static float EasyCurve(float t, float intensity)
    {
        return Mathf.Pow(MathF.Sin(Mathf.Clamp01(t) * MathF.PI / 2f), intensity);
    }

    public static float EasyCurveInverse(float t, float intensity)
    {
        return MathF.Asin(MathF.Pow(Mathf.Clamp01(t), 1f / intensity)) * 2f / MathF.PI;
    }
}

public class SpeedAttribute : PropertyAttribute { }

#if UNITY_EDITOR
[UnityEditor.CustomPropertyDrawer(typeof(SpeedAttribute))]
public class SpeedAttributeDrawer : UnityEditor.PropertyDrawer
{
    public override void OnGUI(Rect position, UnityEditor.SerializedProperty property, GUIContent label)
    {
        const float kphWidth = 100f;

        position.width -= kphWidth;
        UnityEditor.EditorGUI.PropertyField(position, property, label);
        property.floatValue = Mathf.Max(0f, property.floatValue);

        position.x += position.width;
        position.width = kphWidth;
        GUI.Label(position, $"{property.floatValue * 3.6f} Km/h");
    }
}

#endif

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
    public struct Axle
    {
        public WheelData wheelDataL;
        public WheelData wheelDataR;

        public float VisualScale; // TODO: remove, mesh should be properly sized
        public GameObject LeftWheel;
        public GameObject RightWheel;

        public void OnGUI(Camera cam, Transform transform, in CarSettings.Axle axle)
        {
            Vector3 localL = new Vector3(axle.Width * -0.5f, axle.Offset.y, axle.Offset.x);
            Vector3 localR = new Vector3(axle.Width * 0.5f, axle.Offset.y, axle.Offset.x);

            Vector3 wsL = transform.TransformPoint(localL);
            Vector3 wsR = transform.TransformPoint(localR);

            Vector3 screenPos = cam.WorldToScreenPoint(wsL);
            GUI.Label(new Rect(screenPos.x, Screen.height - screenPos.y, 150, 130),
                wheelDataL.compression.ToString("F2"), style);
            screenPos = cam.WorldToScreenPoint(wsR);
            GUI.Label(new Rect(screenPos.x, Screen.height - screenPos.y, 150, 130),
                wheelDataR.compression.ToString("F2"), style);
        }
    }

    public CarSettings Settings;

    [Header("Axles")]
    public Axle Front, Rear;

    [Header("Input")] // TODO: move me to external component
    public bool controllable = true;

    [Header("Debug")]
    public bool debugDraw = true;

    //////////////////////////////////////////////////////////////////////////////////////////////////////

    float steerAngle;
    float steerAngleVelocity;

    float accelerator = 0f; // -1..1, -1 is reverse, 1 is forward
    float brake = 0f; // 0..1

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

        steerAngle = 0f;
    }

    void Start()
    {
        style.normal.textColor = Color.red;
        
        TryGetComponent(out rb);
        rb.centerOfMass = Settings.CenterOfMass;
    }

    void OnValidate()
    {
        //HACK: to apply steering in edit mode
        if (rb == null)
            TryGetComponent(out rb);

        ApplyVisual();
        CalculateAckermannSteering();
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
        float speed = GetSpeed();
        float dt = Time.fixedDeltaTime;

        if (accelerator >= 0f)
            return accelerator * Settings.Forward.GetAccelerationForceMagnitude(speed, dt);
        else
            return accelerator * Settings.Reverse.GetAccelerationForceMagnitude(-speed, dt);
    }

    void Steering(float steeringWheel, float speed)
    {
        speed = Mathf.Abs(speed);
        float steerLimit = Mathf.Lerp(Settings.SteerAngleLimitAtLowSpeed, Settings.SteerAngleLimitAtHighSpeed,
            Utilities.EasyCurve(speed / Settings.MaxSpeed, Settings.SteerAngleLimitCurve));
        float steerTo = steeringWheel * steerLimit;

        steerAngle = Mathf.SmoothDamp(steerAngle, steerTo, ref steerAngleVelocity,
            Settings.SteeringSmooth, float.MaxValue, Time.fixedDeltaTime);
    }

    void UpdateInput()
    {
        float v = 0f;
        float wheel = 0f;
        bool isBrakeNow = false;

        if (controllable)
        {
            v = Input.GetAxis("Vertical");
            wheel = Input.GetAxis("Horizontal");
            if (Input.GetKey(KeyCode.R)) 
                ResetToValidPosition();
            isBrakeNow = Input.GetKey(KeyCode.RightControl) || Input.GetKey(KeyCode.LeftControl);
        }

        brake = 0f;

        float speed = GetSpeed();
        if (speed > Settings.AutoParkThreshold)
        {
            accelerator = Mathf.Max(v, 0f);
            brake = -Mathf.Min(v, 0f);
        }
        else if (speed < -Settings.AutoParkThreshold)
        {
            accelerator = Mathf.Min(v, 0);
            brake = Mathf.Max(v, 0f);
        }
        else
        {
            accelerator = v;
        }

        if (isBrakeNow)
            brake = 1f;

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

        int poweredWheels = 0;
        if (Settings.Front.IsPowered) poweredWheels += 2;
        if (Settings.Rear.IsPowered) poweredWheels += 2;

        CalculateAxleForces(Settings.Front, ref Front, 4, poweredWheels);
        CalculateAxleForces(Settings.Rear, ref Rear, 4, poweredWheels);

        if (TouchingGround())
            Downforce();
        else
            KeepUpwards();
    }

    private bool TouchingGround()
    {
        return Front.wheelDataL.isOnGround  || Front.wheelDataR.isOnGround
            || Rear.wheelDataL.isOnGround || Rear.wheelDataR.isOnGround;
    }

    private void KeepUpwards()
    {
        Vector3 up = transform.worldToLocalMatrix * Vector3.up;

        Vector3 forceUp = new Vector3() { z = up.x * Settings.FlightStabilizationForce * Time.fixedDeltaTime };
        forceUp = transform.TransformDirection(forceUp);
        rb.AddTorque(forceUp, ForceMode.Acceleration);
    }

    private void Downforce()
    {
        Vector3 carDown = transform.TransformDirection(new Vector3(0.0f, -1.0f, 0.0f));

        float speed = Mathf.Abs(GetSpeed());

        float downForceAmount = 1f - Utilities.EasyCurve(1f - speed / Settings.MaxSpeed, Settings.DownForceIntensity);

        rb.AddForce(Settings.DownForce * downForceAmount * carDown, ForceMode.Acceleration);
    }

    void OnGUI()
    {
        if (!controllable)
            return;

        float speed = GetSpeed();
        GUI.Label(new Rect(30, 20, 150, 130), string.Format("{0:F2} km/h", speed * 3.6f));

        GUI.Label(new Rect(30, 50, 1500, 130), $"Steering angle {steerAngle:F2}");
        GUI.Label(new Rect(30, 80, 1500, 130), $"Front Suspensions: {Front.wheelDataL.compression:F2} {Front.wheelDataR.compression:F2}");
        GUI.Label(new Rect(30, 100, 1500, 130), $"Rear Suspensions: {Rear.wheelDataL.compression:F2} {Rear.wheelDataR.compression:F2}");

        Camera cam = Camera.current;
        if (cam == null)
            return;

        if (debugDraw)
        {
            Front.OnGUI(cam, transform, Settings.Front);
            Rear.OnGUI(cam, transform, Settings.Rear);
        }
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

    void CalculateWheelForces(in CarSettings.Axle settings, in Vector3 wsDownDirection, ref WheelData wheelData, in Vector3 wsAttachPoint, int totalWheelsCount, int numberOfPoweredWheels)
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
        float traceLen = settings.RelaxedLength + settings.Radius;

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
            float relaxSpeed = settings.Restitution;
            wheelData.compressionPrev = wheelData.compression;
            wheelData.compression = Mathf.Clamp01(wheelData.compression - dt * relaxSpeed);
            return;
        }

        // Consider wheel radius
        float suspLenNow = wheelData.touchPoint.distance - settings.Radius;

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

        wheelData.compression = 1.0f - Mathf.Clamp01(suspLenNow / settings.RelaxedLength);

        // Hooke's law (springs)
        // F = -k x 

        // Spring force (try to reset compression from spring)
        float springForce = wheelData.compression * -settings.Stiffness;
        suspForceMag += springForce;

        // Damping force (try to reset velocity to 0)
        float suspCompressionVelocity = (wheelData.compression - wheelData.compressionPrev) / dt;
        wheelData.compressionPrev = wheelData.compression;

        float damperForce = -suspCompressionVelocity * settings.Damping;
        suspForceMag += damperForce;

        // Only consider component of force that is along the contact normal.
        float denom = Vector3.Dot(wheelData.touchPoint.normal, -wsDownDirection);
        suspForceMag *= denom;

        // Apply suspension force
        Vector3 suspForce = wsDownDirection * suspForceMag;
        rb.AddForceAtPosition(suspForce, wheelData.touchPoint.point, ForceMode.Acceleration);
        if (debugDraw)
            Debug.DrawRay(wheelData.touchPoint.point, suspForce, Color.yellow);

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
        Vector3 slidingForce = 1f / dt / totalWheelsCount * slideVelocity;

        if (debugDraw)
            Debug.DrawRay(wheelData.touchPoint.point, slideVelocity, Color.red);

        float wheelFriction = Mathf.Clamp01(settings.LateralFriction);
        float groundFriction = wheelData.touchPoint.collider != null
            && wheelData.touchPoint.collider.sharedMaterial != null
            ? wheelData.touchPoint.collider.sharedMaterial.staticFriction : 0.6f;
        float lateralFriction = 0.5f * (wheelFriction + groundFriction);

        // Simulate perfect static friction
        Vector3 frictionForce = slidingForce * -lateralFriction;

        // Remove friction along roll-direction of wheel 
        Vector3 longitudinalForce = Vector3.Dot(frictionForce, c_fwd) * c_fwd;

        bool shouldPark = Mathf.Abs(GetSpeed()) < Settings.AutoParkThreshold;
        bool isAccelerating = Mathf.Abs(accelerationForceMagnitude) > 0.01f;

        // Apply braking force or rolling resistance force or nothing
        if (shouldPark || brake > 0f)
        {
            float mag = longitudinalForce.magnitude;
            float force = mag > 0f ? Mathf.Clamp(settings.BrakeForce, 0.0f, mag) / mag : 0f;
            longitudinalForce *= 1f - brake * force;
        }
        else if (!isAccelerating)
        {
            // Apply rolling-friction (automatic slow-down) only if player don't press the accelerator
            float rollingK = 1.0f - Mathf.Clamp01(settings.RollingFriction);
            longitudinalForce *= rollingK;
        }

        if (debugDraw)
        {
            Debug.DrawRay(wheelData.touchPoint.point, frictionForce, Color.red);
            Debug.DrawRay(wheelData.touchPoint.point, (frictionForce - longitudinalForce), Color.blue);
            Debug.DrawRay(wheelData.touchPoint.point, longitudinalForce, Color.white);
        }

        frictionForce -= longitudinalForce;
            
        // Apply resulting force
        rb.AddForceAtPosition(frictionForce, wheelData.touchPoint.point, ForceMode.Acceleration);

        // Engine force
        if (settings.IsPowered && isAccelerating)
        {
            Vector3 accForcePoint = wheelData.touchPoint.point - (wsDownDirection * 0.2f);
            Vector3 engineForce = c_fwd * accelerationForceMagnitude / numberOfPoweredWheels / dt;
            rb.AddForceAtPosition(engineForce, accForcePoint, ForceMode.Acceleration);

            if (debugDraw)
                Debug.DrawRay(accForcePoint, engineForce, Color.green);
        }
    }

    void CalculateAxleForces(in CarSettings.Axle settings, ref Axle axle, int totalWheelsCount, int numberOfPoweredWheels)
    {
        Vector3 wsDownDirection = transform.TransformDirection(Vector3.down);
        wsDownDirection.Normalize();

        Vector3 localL = new Vector3(settings.Width * -0.5f, settings.Offset.y, settings.Offset.x);
        Vector3 localR = new Vector3(settings.Width * 0.5f, settings.Offset.y, settings.Offset.x);

        Vector3 wsL = transform.TransformPoint(localL);
        Vector3 wsR = transform.TransformPoint(localR);

        CalculateWheelForces(settings, wsDownDirection, ref axle.wheelDataL, wsL, totalWheelsCount, numberOfPoweredWheels);
        CalculateWheelForces(settings, wsDownDirection, ref axle.wheelDataR, wsR, totalWheelsCount, numberOfPoweredWheels);

        // http://projects.edy.es/trac/edy_vehicle-physics/wiki/TheStabilizerBars
        // Apply "stablizer bar" forces
        float travelL = 1.0f - Mathf.Clamp01(axle.wheelDataL.compression);
        float travelR = 1.0f - Mathf.Clamp01(axle.wheelDataR.compression);

        float antiRollForce = (travelL - travelR) * settings.AntiRollForce;
        if (axle.wheelDataL.isOnGround)
        {
            rb.AddForceAtPosition(wsDownDirection * antiRollForce, axle.wheelDataL.touchPoint.point, ForceMode.Acceleration);
            if (debugDraw)
                Debug.DrawRay(axle.wheelDataL.touchPoint.point, wsDownDirection * antiRollForce, Color.green);
        }

        if (axle.wheelDataR.isOnGround)
        {
            rb.AddForceAtPosition(wsDownDirection * -antiRollForce, axle.wheelDataR.touchPoint.point, ForceMode.Acceleration);
            if (debugDraw)
                Debug.DrawRay(axle.wheelDataR.touchPoint.point, wsDownDirection * -antiRollForce, Color.green);
        }
    }

    void CalculateAckermannSteering()
    {
        Debug.Assert(Mathf.Approximately(transform.lossyScale.x, 1f)
            && Mathf.Approximately(transform.lossyScale.y, 1f)
            && Mathf.Approximately(transform.lossyScale.z, 1f));

        // Get turning circle radius for steering angle input
        float turningCircleRadius = Settings.AxleSeparation / Mathf.Tan(steerAngle * Mathf.Deg2Rad);

        // Make front inside tire turn sharper and outside tire less sharp based on turning circle radius
        float steerAngleLeft = Mathf.Atan(Settings.AxleSeparation / (turningCircleRadius + (Settings.FrontWheelsSeparation / 2f)));
        float steerAngleRight = Mathf.Atan(Settings.AxleSeparation / (turningCircleRadius - (Settings.FrontWheelsSeparation / 2f)));

        Front.wheelDataL.yawRad = steerAngleLeft;
        Front.wheelDataR.yawRad = steerAngleRight;
    }

    void CalculateWheelVisualTransform(in Vector3 wsAttachPoint, in Vector3 wsDownDirection, float lengthRelaxed,
        in WheelData data, bool leftWheel, out Vector3 pos, out Quaternion rot)
    {
        float suspCurrentLen = Mathf.Clamp01(1.0f - data.compression) * lengthRelaxed;

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

    void CalculateWheelRotationFromSpeed(CarSettings settings, in CarSettings.Axle axleSettings,
        ref WheelData data, in Vector3 wsPos)
    {
        if (rb == null)
        {
            data.visualRotationRad = 0.0f;
            return;
        }

        float rps;
        if (axleSettings.IsPowered && brake > 0f && Mathf.Abs(accelerator) > 0f)
        {
            rps = accelerator * settings.BurnRotationSpeed;
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
            float wheelLengthMeters = 2 * Mathf.PI * axleSettings.Radius;

            // Wheel "Revolutions per second";
            rps = tireLongSpeed / wheelLengthMeters;
        }
            
        data.rotationsPerSecond = Mathf.SmoothDamp(data.rotationsPerSecond, rps,
            ref data.rotationsPerSecondSpeed, settings.WheelRotationSmoothing);

        float deltaRot = Mathf.PI * 2.0f * data.rotationsPerSecond * Time.deltaTime;

        data.visualRotationRad += deltaRot;
    }

    void ApplyVisual()
    {
        Vector3 wsDownDirection = transform.TransformDirection(Vector3.down);
        wsDownDirection.Normalize();

        ApplyAxleVisual(ref wsDownDirection, ref Front, Settings.Front);
        ApplyAxleVisual(ref wsDownDirection, ref Rear, Settings.Rear);
    }

    private void ApplyAxleVisual(ref Vector3 wsDownDirection, ref Axle axle, in CarSettings.Axle axleSettings)
    {
        Vector3 localL = new Vector3(axleSettings.Width * -0.5f, axleSettings.Offset.y, axleSettings.Offset.x);
        Vector3 localR = new Vector3(axleSettings.Width * 0.5f, axleSettings.Offset.y, axleSettings.Offset.x);

        Vector3 wsL = transform.TransformPoint(localL);
        Vector3 wsR = transform.TransformPoint(localR);

        Vector3 wsPos;
        Quaternion wsRot;

        Vector3 wheelScale = axleSettings.Radius * axle.VisualScale * Vector3.one; // TODO: remove

        if (axle.LeftWheel != null)
        {
            CalculateWheelVisualTransform(wsL, wsDownDirection, axleSettings.RelaxedLength, axle.wheelDataL,
                true, out wsPos, out wsRot);
            axle.LeftWheel.transform.SetPositionAndRotation(wsPos, wsRot);
            axle.LeftWheel.transform.localScale = wheelScale;

            CalculateWheelRotationFromSpeed(Settings, axleSettings, ref axle.wheelDataL, wsPos);
        }

        if (axle.RightWheel != null)
        {
            CalculateWheelVisualTransform(wsR, wsDownDirection, axleSettings.RelaxedLength, axle.wheelDataR,
                false, out wsPos, out wsRot);
            axle.RightWheel.transform.SetPositionAndRotation(wsPos, wsRot);
            axle.RightWheel.transform.localScale = wheelScale;

            CalculateWheelRotationFromSpeed(Settings, axleSettings, ref axle.wheelDataR, wsPos);
        }
    }

#if UNITY_EDITOR
    void OnDrawGizmosAxle(Vector3 wsDownDirection, in CarSettings.Axle axleSettings, in Axle axle)
    {
        Vector3 localL = new Vector3(axleSettings.Width * -0.5f, axleSettings.Offset.y, axleSettings.Offset.x);
        Vector3 localR = new Vector3(axleSettings.Width * 0.5f, axleSettings.Offset.y, axleSettings.Offset.x);

        Vector3 wsL = transform.TransformPoint(localL);
        Vector3 wsR = transform.TransformPoint(localR);

        Gizmos.color = Color.magenta;

        //draw axle
        Gizmos.DrawLine(wsL, wsR);

        //draw line to com
        Gizmos.DrawLine(transform.TransformPoint(
            new Vector3(0.0f, axleSettings.Offset.y, axleSettings.Offset.x)),
            transform.TransformPoint(Settings.CenterOfMass));

        DrawWheelGizmo(ref wsDownDirection, axleSettings, true, axle.wheelDataL, wsL);
        DrawWheelGizmo(ref wsDownDirection, axleSettings, false, axle.wheelDataR, wsR);
    }

    private void DrawWheelGizmo(ref Vector3 wsDownDirection, in CarSettings.Axle axle,
        bool leftWheel, in WheelData wheelData, in Vector3 wsFrom)
    {
        Gizmos.color = wheelData.isOnGround ? Color.yellow : Color.red;
        UnityEditor.Handles.color = Gizmos.color;

        float suspCurrentLen = Mathf.Clamp01(1.0f - wheelData.compression) * axle.RelaxedLength;

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
        Gizmos.DrawLine(wsTo, wsTo + wsForward * axle.Radius);

        // Draw wheel
        UnityEditor.Handles.DrawWireDisc(wsTo, wsAxle, axle.Radius);
        UnityEditor.Handles.DrawWireDisc(wsTo + wsAxle * wheelWidth, wsAxle, axle.Radius);
        UnityEditor.Handles.DrawWireDisc(wsTo - wsAxle * wheelWidth, wsAxle, axle.Radius);
    }

    void OnDrawGizmos()
    {
        Vector3 wsDownDirection = transform.TransformDirection(Vector3.down);
        wsDownDirection.Normalize();
        OnDrawGizmosAxle(wsDownDirection, Settings.Front, Front);
        OnDrawGizmosAxle(wsDownDirection, Settings.Rear, Rear);
        Gizmos.DrawSphere(transform.TransformPoint(Settings.CenterOfMass), 0.1f);
    }
#endif
}

