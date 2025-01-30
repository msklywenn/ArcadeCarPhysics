using System;
using UnityEngine;

[CreateAssetMenu]
public class CarSettings : ScriptableObject
{
    [Serializable]
    public struct Axle
    {
        [Header("Axle settings")]
        [Min(0f)] public float Width;
        public Vector2 Offset;
        public bool IsPowered;
        
        [Header("Wheel settings")]
        [Min(0.01f)] public float Radius;

        [Tooltip("Tire laterial friction normalized to 0..1")]
        [Range(0.0f, 1.0f)] public float LateralFriction; // TODO use static friction from physics material

        [Range(0.0f, 0.1f)] public float RollingFriction;

        [Tooltip("Brake force magnitude")]
        [Min(0.01f)] public float BrakeForce;

        //--
        [Header("Suspension settings")]
        [Tooltip("Suspension Stiffness (Suspension 'Power'")]
        [Min(0.01f)] public float Stiffness;

        [Tooltip("Suspension Damping (Suspension 'Bounce')")]
        [Min(0.01f)] public float Damping;

        [Tooltip("Suspension Restitution")]
        [Min(0.01f)] public float Restitution;

        [Min(0.01f)] public float RelaxedLength;

        [Tooltip("Stabilizer bar anti-roll force")]
        [Min(0.01f)] public float AntiRollForce;
    }

    [Serializable]
    public struct Engine
    {
        [Speed] public float TopSpeed;

        [Tooltip("How quickly to reach top speed. <1 is quickly, 1 is linearly, >1 is slowly")]
        [Min(0.01f)] public float Acceleration;

        [Tooltip("Time to reach top speed")]
        [Min(0.01f)] public float Duration;

        public float GetAccelerationForceMagnitude(float speed, float dt)
        {
            float t = Duration * Utilities.EasyCurveInverse(speed / TopSpeed, Acceleration);
            float desiredSpeed = TopSpeed * Utilities.EasyCurve((t + dt) / Duration, Acceleration);
            desiredSpeed = Mathf.Clamp(desiredSpeed, 0f, TopSpeed);
            float acceleration = desiredSpeed - speed;

            float intensity = acceleration;
            intensity = Mathf.Max(intensity, 0.0f);
            return intensity;
        }
    };

    public Vector3 CenterOfMass = new Vector3(0f, -0.48f, 0.38f);

    [Tooltip("Only used as a reference for downforce and steering settings")]
    public float MaxSpeed = 200f; // used only for settings

    public Engine Forward = new Engine { TopSpeed = 130 / 3.6f, Acceleration = 0.8f, Duration = 10f };
    public Engine Reverse = new Engine { TopSpeed = 40 / 3.6f, Acceleration = 0.8f, Duration = 3f };

    [Header("Steering")]
    public float SteeringSpeed = 60f;
    public float SteeringResetSpeed = 50f;
    // x - speed in km/h
    // y - angle in degrees
    public float SlowSteerAngleLimit = 40f;
    public float FastSteerAngleLimit = 4f;
    public float SlowFastSteerCurve = 0.9f;//TODO: a better name...

    [Header("Other")]
    [Tooltip("Stabilization in flight (torque)")]
    public float FlightStabilizationForce = 8.0f;

    [Tooltip("Stabilization in flight (Ang velocity damping)")]
    public float FlightStabilizationDamping = 0.0f;

    public float DownForceIntensity = 0.7f;

    [Tooltip("Downforce")]
    public float DownForce = 10f;

    [Tooltip("Rotation speed when braking and accelerating at the same time")]
    public float BurnRotationSpeed = 5f;
    public float WheelRotationSmoothing = 0.1f;

    public Axle Front = new Axle
    {
        Width = 1.55f,
        Offset = new Vector2(1.51f, 0f),
        IsPowered = false,
        Radius = 0.3f,
        LateralFriction = 0.6f,
        RollingFriction = 0.005f,
        BrakeForce = 14,
        Stiffness = 10000,
        Damping = 3000,
        Restitution = 1f,
        RelaxedLength = 0.45f,
        AntiRollForce = 10000,
    };
    public Axle Rear = new Axle
    {
        Width = 1.55f,
        Offset = new Vector2(-1.29f, 0f),
        IsPowered = true,
        Radius = 0.3f,
        LateralFriction = 0.6f,
        RollingFriction = 0.005f,
        BrakeForce = 2,
        Stiffness = 8000,
        Damping = 3000,
        Restitution = 1f,
        RelaxedLength = 0.45f,
        AntiRollForce = 10000,
    };

    [SerializeField] public float AxleSeparation;
    [SerializeField] public float FrontWheelsSeparation;

    void OnValidate()
    {
        // Calculate our chassis assuming uniform scale
        Vector3 frontOffset = new Vector3(0.0f, Front.Offset.y, Front.Offset.x);
        Vector3 rearOffset = new Vector3(0.0f, Rear.Offset.y, Rear.Offset.x);
        AxleSeparation = Vector3.Distance(frontOffset, rearOffset);

        Vector3 leftWheel = new Vector3(Front.Width * -0.5f, Front.Offset.y, Front.Offset.x);
        Vector3 rightWheel = new Vector3(Front.Width * 0.5f, Front.Offset.y, Front.Offset.x);
        FrontWheelsSeparation = Vector3.Distance(leftWheel, rightWheel);

    }
}

