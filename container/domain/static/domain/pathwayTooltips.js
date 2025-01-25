function makePathwayTooltip(pathwayName) {
    return `
    <div class="pathway-tooltip">
      <p style="text-align: center">Click to inspect pathway: ${pathwayName}</p>
      <button class="btn btn-light m-1">Analyse</button>
      <button class="btn btn-light m-1">Visualise</button>
    </div>
    `;
}

function addTooltip(triggers) {
    triggers.forEach(function (trigger) {
        // get pathway name from the trigger title
        const pathwayName = trigger.title;
        // add the tooltip
        trigger.innerHTML += makePathwayTooltip(pathwayName);

        // positioning the tooltip
        trigger.addEventListener('mouseover', function() {
            const tooltip = document.querySelector('.pathway-tooltip');
            const parentRect = this.getBoundingClientRect();
            tooltip.style.display = 'hidden';
            tooltip.style.left = `${parentRect.left + parentRect.width / 2}px`;
            tooltip.style.top = `${parentRect.top - tooltip.offsetHeight - 10}px`;
            console.log(parentRect.left, parentRect.width, parentRect.top, tooltip.offsetHeight);
            tooltip.style.display = 'block';
        });

        trigger.addEventListener('mouseout', function() {
            const tooltip = document.querySelector('.pathway-tooltip');
            setTimeout(function () {
                tooltip.style.display = 'none';
            }, 100)
        });
    });
}

