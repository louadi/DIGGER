function makePathwayTooltip(pathwayName) {
    return `
    <div class="pathway-tooltip">
      <p style="text-align: center">Click to inspect pathway: ${pathwayName}</p>
      <button class="btn btn-light m-1">Analyse</button>
      <button class="btn btn-light m-1">Visualise</button>
    </div>
    `;
}

function tooltipTransitions(triggers, pathwayName) {
    triggers.forEach(function (trigger) {
        // add the tooltip
        trigger.innerHTML += makePathwayTooltip(pathwayName);

        // add the hover events
        trigger.addEventListener('mouseenter', function () {
            var tooltip = trigger.querySelector('.pathway-tooltip');
            tooltip.style.display = 'block';
            setTimeout(function () {
                tooltip.classList.add('hovered');
            }, 10); // Small delay to trigger transition
        });

        trigger.addEventListener('mouseleave', function () {
            var tooltip = trigger.querySelector('.pathway-tooltip');
            setTimeout(function () {
                if (!tooltip.matches(':hover')) {
                    tooltip.classList.remove('hovered');
                    setTimeout(function () {
                        tooltip.style.display = 'none';
                    }, 100); // Match the transition duration
                }
            }, 100); // Delay for smooth fade-out effect
        });

        var tooltip = trigger.querySelector('.pathway-tooltip');
        tooltip.addEventListener('mouseenter', function () {
            tooltip.classList.add('hovered');
        });

        tooltip.addEventListener('mouseleave', function () {
            setTimeout(function () {
                if (!trigger.matches(':hover')) {
                    tooltip.classList.remove('hovered');
                    setTimeout(function () {
                        tooltip.style.display = 'none';
                    }, 100); // Match the transition duration
                }
            }, 100); // Delay for smooth fade-out effect
        });
    });
}