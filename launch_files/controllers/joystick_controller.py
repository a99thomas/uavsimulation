import pygame


class Joystick:
    def __init__(self):
        pygame.init()
        pygame.joystick.init()

        if pygame.joystick.get_count() == 0:
            print("No joystick detected")
            self.joystick = None
        else:
            self.joystick = pygame.joystick.Joystick(0)
            self.joystick.init()

    def get_transmitter_values(self):
        pygame.event.pump()
        if self.joystick:
            throttle = self.remap_value(self.joystick.get_axis(2), -0.71, 0.71, 0, 1)
            rudder = self.remap_value(self.joystick.get_axis(3), -0.71, 0.71, 0.25, -0.25)
            aileron = self.remap_value(self.joystick.get_axis(0), -0.71, 0.71, -0.25, 0.25)
            elevator = self.remap_value(self.joystick.get_axis(1), -0.71, 0.71, 0.25, -.4)
            trigger = self.joystick.get_button(1)
            if trigger == 1:
                trigger = 0
            else:
               trigger = 1
            return throttle, rudder, aileron, elevator, trigger
        else:
            return None, None, None, None, None

    def remap_value(self, value, old_min, old_max, new_min, new_max):
        # Handle cases where value is outside old range
        if value < old_min:
            value = old_min
        elif value > old_max:
            value = old_max
        old_range = old_max - old_min
        new_range = new_max - new_min
        scaled_value = (value - old_min) / old_range
        new_value = scaled_value * new_range + new_min
        return new_value


def main():
  joystick = Joystick()

  running = True
  while running:
    for event in pygame.event.get():
      if event.type == pygame.QUIT:
        running = False

    throttle, rudder, aileron, elevator, trigger = joystick.get_transmitter_values()
    # print(trigger)
    # if throttle is not None:
    #     # ... your code here to process joystick values ...

    # else:
    #     print("Failed to initialize joystick")
    #     break

  pygame.quit()

if __name__ == "__main__":
  main()
